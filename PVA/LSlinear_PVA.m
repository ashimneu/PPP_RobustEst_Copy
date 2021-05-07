function [pos,msr_res,GDOP,nsv,dnsv,nprior,dnprior,x_post,P_post,delta_x,Jydiag,by,H_pos_vel,res_std] = ...
    LSlinear_PVA(p,cpt,grdpos,x_prior,P_prior_in)

 
% ind_s = find(cpt.num_sv([1,3,4]) ~= 0); %<---
ind_s = [1 2 3];
xk = [x_prior(1:9);x_prior(9+ind_s);x_prior(13)]; % point of linearization
ind_no = find(cpt.num_sv([1,3,4]) == 0);
Pcov = P_prior_in;
% Pcov(9+ind_no,:) = []; %<---
% Pcov(:,9+ind_no) = []; %<---

y_R = cpt.corr_range;
y_D = cpt.dp_range;
num = length(y_R); % number of measurement
H = zeros(num,3);
% num_clk = sum(cpt.num_sv~=0); %<---
% H_clk = zeros(num,num_clk); %<---
% ind = find(cpt.num_sv~=0); %<---
% start = 1;
% for i = 1:num_clk
%     H_clk(start:start+cpt.num_sv(ind(i))-1,i)=1;
%     start = start + cpt.num_sv(ind(i));
% end
%--------------------------------------------------------------------------
H_clk = zeros(num,3);
ind = find(cpt.num_sv([1 3 4])~=0); % indices of visibile Constellations
num_sv = cpt.num_sv(find(cpt.num_sv~=0)); % # SV of each visibile Constellations
start = 1; 
for constellation = 1:1:numel(ind)
    H_clk(start:start+num_sv(constellation)-1,ind(constellation))=1;
    start = start + num_sv(constellation);
end
%--------------------------------------------------------------------------
H_clk(:,1) = 1;
R = zeros(num,1);  
r = zeros(num,1);
rv = zeros(num,1);
clk_bia = zeros(num,1);
lamda = p.c*[ones(cpt.num_sv(1),1)/p.L1freq; ones(cpt.num_sv(2),1)/p.L1freq; ...
    ones(cpt.num_sv(3),1)/p.E1freq; ones(cpt.num_sv(4),1)/p.B1freq];
y_D = -lamda.*y_D;
if p.post_mode == 1
    if p.IGS_enable == 1
        s_pos_ecef = cpt.s_pos_prc;
        s_v_ecef = cpt.s_v_ecef;
    else
        s_pos_ecef = cpt.s_pos_ecef;
        s_v_ecef = cpt.s_v_ecef;
    end
else
    s_pos_ecef = cpt.s_pos_ecef;
    s_v_ecef = cpt.s_v_ecef;
end
for iter=1:1
for j=1:num
    R(j)=norm(s_pos_ecef(:,j)-xk(1:3));
    H(j,:) = (xk(1:3)-s_pos_ecef(:,j))'/R(j)+...
        [-s_pos_ecef(2,j)*p.omge/p.c s_pos_ecef(1,j)*p.omge/p.c 0];
    rv(j) = -H(j,:)*(s_v_ecef(:,j)-xk(4:6))+...
        sagnac_v(p,[s_pos_ecef(:,j);s_v_ecef(:,j)],xk);
    r(j) = R(j)+sagnac(p,s_pos_ecef(:,j),xk);
    ind = find(H_clk(j,:)==1);
    clk_bia(j) = sum(xk(9+ind));
end
H_os = [H,zeros(size(H)),zeros(num,3),H_clk,zeros(num,1);
    zeros(size(H)),H,zeros(num,3),zeros(size(H_clk)),ones(num,1)];
res_R = y_R - r - clk_bia + cpt.outliervec;
res_D = (y_D - (rv+xk(end))) + cpt.outliervec;
yCov = blkdiag(p.sig_y^2.*eye(num),p.sig_y_dop^2.*eye(num)); % noise covariance
resCov = H_os*P_prior_in*H_os' + yCov; % Section VI:A, 1st sentence.
res_std = sqrt(diag(resCov));

%  Innovation Covariance
S = H_os*Pcov*H_os' + yCov;
% Optimal Kalman Gain
Kk = Pcov*H_os'*S^(-1);
% Posterior state estimate
dx = Kk*[res_R;res_D];
xk = xk + dx;
% Posterior state covariance estimate
T = Kk*H_os;
Pcov = (eye(size(T))-T)*Pcov;
end
% p.state_PVA([1:9,9+ind_s,13]) = xk;
% ivec = [1:9,9+ind_s,13];
% for l = 1:length(ivec)
%     p.P_cov(ivec(l),[1:9,9+ind_s,13]) = Pcov(l,:);
% end
% end

%------------------------%
pos = xk(1:3);
x_post = xk; % Posterior xhat
err1 = norm(grdpos - pos); % norm(p.P_base - pos); % error check for debugging

%--------------------------------------------------------------------------
ind_s2 = find(cpt.num_sv([1,3,4]) ~= 0);
H_clk2 = H_clk(:,ind_s2);
H_os2 = [H,zeros(size(H)),H_clk2,zeros(num,1);
            zeros(size(H)),H,zeros(size(H_clk2)),ones(num,1)];
delta_x2 = ((H_os2'*H_os2)^(-1)) * (H_os2'*[res_R;res_D]);
%------------------------%
x_prior2 = [x_prior(1:6); x_prior(9+ind_s2); x_prior(end)];
x_hat2 = x_prior2 + delta_x2;
pos2 = x_prior(1:3) + delta_x2(1:3);
err2 = norm(grdpos - pos2);


% Posterior xhat covariance
iR_psr = eye(num)./p.sig_y^2; 
iR_dop = eye(num)./p.sig_y_dop^2;
invR   = blkdiag(iR_psr,iR_dop);
P_post = (H_os'*invR*H_os + P_prior_in^(-1))^(-1);

% Posterior information (obtained from observations only)
Jydiag = diag(H_os'*yCov^(-1)*H_os);
 
%--------------------------%
% Measurement Residual
delta_x = dx;
mu_x = zeros(size(x_prior));
E_R = sqrtm(yCov^(-1));
E_P = sqrtm(P_prior_in^(-1));
A = [E_R*H_os; E_P];
c = [E_R*[res_R;res_D]; E_P*mu_x];
residual = c - A*delta_x;
msr_res  = residual(1:2*num);

%--------------------------%
% Compute GDOP
GDOP = sqrt(trace((H'*iR_psr*H)^(-1)));

%-----------------------%
n       = size(x_prior);
nsv     = 2*num; % number of measurements used
dnsv    = 0; % number of measurements discarded
nprior  = size(n,1); % number of priors used
dnprior = 0; % number of priors discarded
by      = ones(2*num,1);
H_pos_vel = H_os(:,1:6);

% %---------------------------%
% % Compute pseudorange & doppler residual to compute Measurement Covar. R
% ind = find(cpt.num_sv([1 3 4])~=0);;
% linearization_point = [grdpos; zeros(3,1); zeros(3,1); x_post(10+ind); 142.1628]; % true state vector
% for j=1:num
%     R(j)=norm(s_pos_ecef(:,j)-linearization_point(1:3));
%     H(j,:) = (linearization_point(1:3)-s_pos_ecef(:,j))'/R(j)+...
%         [-s_pos_ecef(2,j)*p.omge/p.c s_pos_ecef(1,j)*p.omge/p.c 0];
%     rv(j) = -H(j,:)*(s_v_ecef(:,j)-linearization_point(4:6))+...
%         sagnac_v(p,[s_pos_ecef(:,j);s_v_ecef(:,j)],linearization_point);
%     r(j) = R(j)+sagnac(p,s_pos_ecef(:,j),linearization_point);
%     ind = find(H_clk(j,:)==1);
% %     clk_bia_tru(j) = sum(linearization_point(9+ind));
% end
% % H_os = [H,zeros(size(H)),zeros(num,3),H_clk,zeros(num,1);
% %     zeros(size(H)),H,zeros(num,3),zeros(size(H_clk)),ones(num,1)];
% res_R_tru = y_R - r - clk_bia;
% res_D_tru = (y_D - (rv+linearization_point(end)));
% msr_res_tru = [res_R_tru; res_D_tru];

end