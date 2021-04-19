function [pos,msr_res,GDOP,nsv,dnsv,nprior,dnprior,x_post,P_post,dx,Jydiag,by,H_pos_vel] = ...
    LTSlinear_PVA_useallprior(p,cpt,grdpos,x_prior,P_prior,Option)
ind_s = [1 2 3];
xk = [x_prior(1:9);x_prior(9+ind_s);x_prior(13)]; % point of linearization
ind_no = find(cpt.num_sv([1,3,4]) == 0);
Pcov = P_prior;

y_R = cpt.corr_range;
y_D = cpt.dp_range;
num = length(y_R); % number of measurement
%--------------------------------------------------------------------------
H_clk = zeros(num,3);
ind = find(cpt.num_sv([1 3 4])~=0); % indices of visibile Constellations
num_sv = cpt.num_sv(find(cpt.num_sv~=0)); % # SV of each visibile Constellations
start = 1; 
for constellation = 1:1:numel(ind)
    H_clk(start:start+num_sv(constellation)-1,ind(constellation))=1;
    start = start + num_sv(constellation);
end
H_clk(:,1) = 1;
H = zeros(num,3);
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
res_R = y_R - r - clk_bia;
res_D = (y_D - (rv+xk(end)));

yCov = blkdiag(p.sig_y^2.*eye(num),p.sig_y_dop^2.*eye(num)); % noise covariance
E_R  = chol(yCov^(-1));
% E_P  = chol(Pcov^(-1));
n    = numel(x_prior);
mu_x = zeros(n,1);
A   = E_R*H_os;            %A2 = E_P;      A = [A1; A2];
c   = E_R*[res_R;res_D] ;  %c2 = E_P*mu_x; c = [c1 ; c2];

A1 = [A(1:num,1:3) A(1:num,10:end-1)];     c1 = c(1:num);
A2 = [A(num+1:end,4:6) A(num+1:end,end)];  c2 = c(num+1:end);
A1bar   = A1;
zerocol = ~any(A1,1);  % assigns 0 if col of 0s is present
A1bar(:,zerocol) = []; % removes cols of 0s
A2bar   = A2;
zerocol = ~any(A2,1);  % assigns 0 if col of 0s is present
A2bar(:,zerocol) = []; % removes cols of 0s

switch Option
case 1
    Abar = A;
    zerocol = ~any(A,1);  % assigns 0 if col of 0s is present
    Abar(:,zerocol) = []; % removes cols of 0s
    [rew,~] = ltsregres(Abar,c,'plots',0,'intercept',0);
    by = rew.flag;
case 2
    Abar = A;
    zerocol = ~any(A,1);  % assigns 0 if col of 0s is present
    Abar(:,zerocol) = []; % removes cols of 0s
    rew = lts(Abar,c);
    by = rew.flag;
case 3
%     h = floor((num+n+1)/2)+1; % for highest breakdown point
    % [rew,~] = ltsregres(A1,c1,'plots',0,'intercept',0,'h',h);
    rew1 = lts(A1bar,c1);
    rew2 = lts(A2bar,c2);
    by1  = rew1.flag;    % measurement selection binary vector
    by2  = rew2.flag;    % measurement selection binary vector
    by = [by1; by2];
case 4
    rew1 = lts(A1bar,c1);
%     rew2 = lts(A2bar,c2);
    by1  = rew1.flag;    % measurement selection binary vector
%     by2  = rew2.flag;    % measurement selection binary vector
    by = [by1; by1];
case 5
%     rew1 = lts(A1bar,c1);
    rew2 = lts(A2bar,c2);
%     by1  = rew1.flag;    % measurement selection binary vector
    by2  = rew2.flag;    % measurement selection binary vector
    by = [by2; by2];
case 6
    h = size(A1bar,1)-1;
    rew = lts(A1bar,c,h);
end 

% by = ones(2*num,1);
bx      = ones(n,1);   % prior selection binary vector; we use all priors
Pby  = diag(by);
Pbx  = eye(n);
PhiH = Pby*H_os;
H_pos_vel = PhiH(:,1:6);

dx = ((PhiH'*yCov^(-1)*PhiH+ P_prior^(-1))^(-1))*...
    (PhiH'*yCov^(-1)*Pby*[res_R;res_D] + P_prior^(-1)*mu_x);
%-----------------------%
xk  = xk + dx;
x_post = xk;
pos = xk(1:3);
err1 = norm(grdpos - pos);
end

%--------------------------------------------------------------------------
ind_s2 = find(cpt.num_sv([1,3,4]) ~= 0);
H_clk2 = H_clk(:,ind_s2);
H_os2 = [H,zeros(size(H)),H_clk2,zeros(num,1);
            zeros(size(H)),H,zeros(size(H_clk2)),ones(num,1)];
delta_x2 = ((H_os2'*H_os2)^(-1)) * (H_os2'*[res_R;res_D]);
%------------------------%
x_prior2 = [x_prior(1:6); x_prior(9+ind_s2); x_prior(end)];
pos2 = x_prior(1:3) + delta_x2(1:3);
err2 = norm(grdpos - pos2);


% Compute posterior xhat covariance
Hbar = PhiH; % Hbar = [PhiH(:,1:10), PhiH(:,end)];
iR_psr = eye(num)./p.sig_y^2; 
iR_dop = eye(num)./p.sig_y_dop^2;
invR   = blkdiag(iR_psr,iR_dop);
P_post = (Hbar'*invR*Hbar + P_prior^(-1))^(-1);

%--------------------------%
% Compute posterior information (obtained from observations only)
Jydiag = diag(PhiH'*invR*PhiH);
 
%--------------------------%
% Compute residual
mu_x = zeros(size(x_prior));
E_R = sqrtm(yCov^(-1));
E_P = sqrtm(P_prior^(-1));
Ab = [E_R*Pby*H_os; E_P*Pbx];
cb = [E_R*Pby*[res_R;res_D]; E_P*Pbx*mu_x];
residual = cb - Ab*dx ;
msr_res  = residual(1:2*num);

% Compute GDOP
Hbar_psr = remove0colrow(PhiH(1:num,1:3),"column");
GDOP = sqrt(trace((Hbar_psr'*iR_psr*Hbar_psr)^(-1)));

nsv     = sum(by);     % count of measurements used
dnsv    = 2*num - nsv; % count of measurements discarded
nprior  = sum(bx);     % count of priors used
dnprior = n - nprior;  % count of priors discarded

end