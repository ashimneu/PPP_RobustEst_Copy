function [pos,msr_res,GDOP,nsv,dnsv,nprior,dnprior,x_post,P_post,delta_x,Jydiag,by,H_pos_vel] = ...
    LSlinear_PVA(p,cpt,grdpos,x_prior,P_prior_in)

ind_s = find(cpt.num_sv([1,3,4]) ~= 0);
xk = [x_prior(1:9);x_prior(9+ind_s);x_prior(13)]; % point of linearization
ind_no = find(cpt.num_sv([1,3,4]) == 0);
Pcov = P_prior_in;
Pcov(9+ind_no,:) = [];
Pcov(:,9+ind_no) = [];

y_R = cpt.corr_range;
y_D = cpt.dp_range;
num = length(y_R); % number of measurement
num_clk = sum(cpt.num_sv~=0);
H = zeros(num,3);
H_clk = zeros(num,num_clk);
ind = find(cpt.num_sv~=0);
start = 1;
for i = 1:num_clk
    H_clk(start:start+cpt.num_sv(ind(i))-1,i)=1;
    start = start + cpt.num_sv(ind(i));
end
H_clk(:,1) = 1;
R = zeros(num,1);  
r = zeros(num,1);
rv = zeros(num,1);
clk_bia = zeros(num,1);
lamda = p.c*[ones(cpt.num_sv(1),1)/p.L1freq;ones(cpt.num_sv(2),1)/p.L1freq;...
    ones(cpt.num_sv(3),1)/p.E1freq;ones(cpt.num_sv(4),1)/p.B1freq];
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
res_R = y_R - r - clk_bia;
res_D = (y_D - (rv+xk(end)));
yCov = eye(num*2); % p.sig_y^2.*eye(2*num); % noise covariance
%  Innovation Covariance
S = H_os*Pcov*H_os' + yCov; %+eye(num*2);
% Optimal Kalman Gain
Kk = Pcov*H_os'*S^(-1);
% Posterior state estimate
dx = Kk*[res_R;res_D];
xk = xk + dx;
% Posterior state covariance estimate
T = Kk*H_os;
Pcov = (eye(size(T))-T)*Pcov;
% end
% p.state_PVA([1:9,9+ind_s,13]) = xk;
% ivec = [1:9,9+ind_s,13];
% for l = 1:length(ivec)
%     p.P_cov(ivec(l),[1:9,9+ind_s,13]) = Pcov(l,:);
% end
% end

H_pos_vel = H_os(:,1:6);
 
% %-----------------------%      
% for j=1:m
%     Range(j)=norm(s_pos_ecef(:,j)-xk(1:3));
%     V= (xk(1:3)-s_pos_ecef(:,j))'/Range(j);
%     H(j,:)=[V 1];  
%     r(j) = Range(j)+sagnac(p,s_pos_ecef(:,j),xk);
%     if ~isempty(H_offset)
%         ind = find(H_offset(j,:)==1);
%         if ~isempty(ind)
%             off(j) = xk(4+ind);
%         end
%     end
% end
% res = y - r - xk(4)-off; % pseudorange residual
% 
% H_pos   = H(:,1:3);   % H for position
% H_vel   = zeros(m,3); % H for velocity
% H_accl  = zeros(m,3); % H for acceleration
% H_c     = [H(:,4), H_offset, zeros(m,1)]; % H for clock bias & drift
% H_os    = [H_pos, H_vel, H_accl, H_c];
% [x_prior,P_prior] = adjust_entries(p,x_prior,P_prior_in,0);

%-----------------------%
n       = size(x_prior);
nsv     = 2*num; % number of measurements used
dnsv    = 0; % number of measurements discarded
nprior  = size(n,1); % number of priors used
dnprior = 0; % number of priors discarded
by      = ones(2*num,1);

yCov = eye(2*num); % p.sig_y^2.*eye(num); % noise covariance
% delta_x = (H_os'*R^(-1)*H_os+ P_prior^(-1))^(-1)*(H_os'*R^(-1)*res);
delta_x = dx;
%------------------------%
% x_hat = x_prior + delta_x;
x_hat = xk;
pos = x_hat(1:3);
clock_bias = x_hat(4);
x_post = x_hat; % Posterior xhat
err1 = norm(grdpos - pos); % norm(p.P_base - pos); % error check for debugging

H_os2 = [H,zeros(size(H)),H_clk,zeros(num,1);
            zeros(size(H)),H,zeros(size(H_clk)),ones(num,1)];
delta_x2 = ((H_os2'*H_os2)^(-1)) * (H_os2'*[res_R;res_D]);
%------------------------%
x_prior2 = [x_prior(1:6); x_prior(9+ind_s); x_prior(end)];
x_hat2 = x_prior2 + delta_x2;
pos2 = x_prior(1:3) + delta_x2(1:3);
err2 = norm(grdpos - pos2);


% Posterior xhat covariance
yCov   = eye(2*num).*p.sig_y^2; % R inverse
P_post = (H_os'*yCov^(-1)*H_os + P_prior_in^(-1))^(-1);

% Posterior information (obtained from observations only)
Jydiag = diag(H_os'*yCov^(-1)*H_os);
 
%--------------------------%
% Measurement Residual
mu_x = zeros(size(x_prior));
E_R = sqrtm(yCov^(-1));
E_P = sqrtm(P_prior_in^(-1));
A = [E_R*H_os; E_P];
c = [E_R*[res_R;res_D]; E_P*mu_x];
residual = c - A*delta_x;
msr_res  = residual(1:2*num);

%--------------------------%
% Compute GDOP
H_os3 = [H,zeros(size(H)); zeros(size(H)),H]; %[H, H_clk];
yCov2 = p.sig_y^2.*eye(2*num); % noise covariance
% delta_x = ((H_os2'*R2^(-1)*H_os2+Pcov^(-1))^(-1))*(H_os2'*R2^(-1)*res);
GDOP = sqrt(trace((H_os3'*yCov2^(-1)*H_os3)^(-1)));

end