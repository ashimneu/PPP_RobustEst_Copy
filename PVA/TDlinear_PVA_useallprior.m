function [pos,msr_res,GDOP,nsv,dnsv,nprior,dnprior,x_post,P_post,dx,Jydiag,by,H_pos_vel] = ...
    TDlinear_PVA_useallprior(p,cpt,grdpos,x_prior,P_prior,TDLambda)

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

%-----------------------%
n       = numel(xk);
mu_x    = zeros(n,1);
yCov = blkdiag(p.sig_y^2.*eye(num),p.sig_y_dop^2.*eye(num)); % noise covariance
R       = p.sig_y^2.*eye(num);
reshat  = H_os*mu_x;
resdiff = abs([res_R; res_D] - reshat); %residual

resCov = H_os*P_prior*H_os' + yCov; % Section VI:A, 1st sentence.
sig_ri = sqrt(diag(resCov));
ratio  = resdiff./sig_ri;
Lvec   = TDLambda.*ones(2*num,1); % vector of thresholds
by     = zeros(2*num,1); % binary decision vector

for i = 1:1:num
    % generate selection vector b_TD
    if ratio(i) <= Lvec(i) 
        by(i) = 1;     % keep res(i)
    else
        by(i) = 0;     % flag res(i)
    end
end

bx   = ones(n,1);   % prior selection binary vector; use all priors
nsv  = sum(by);     % number of measurements used
dnsv = 2*num - nsv; % number of measurements discarded
nprior  = n;        % number of priors used
dnprior = 0;        % number of priors discarded

Pby  = diag(by);
Pbx  = eye(n);
PhiH = Pby*H_os;
H_pos_vel = PhiH(:,1:6);

dx   = ((PhiH'*yCov^(-1)*PhiH+ P_prior^(-1))^(-1))*...
    (PhiH'*yCov^(-1)*Pby*[res_R;res_D] + P_prior^(-1)*mu_x);
%-----------------------%
xk = xk + dx;
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
% Rbar = eye(size(Hbar,1))./p.sig_y^2;
Rbar_R = eye(num)./p.sig_y^2; 
Rbar_D = eye(num)./p.sig_y_dop^2;
Rbar   = blkdiag(Rbar_R,Rbar_D);
P_post = (Hbar'*Rbar*Hbar + P_prior^(-1))^(-1);

%--------------------------%
% Compute posterior information (obtained from observations only)
Jydiag = diag(PhiH'*Rbar*PhiH);

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
PhiH2 = PhiH(:,1:3);
yCov2 = p.sig_y^2.*eye(2*num); % noise covariance
GDOP = sqrt(trace((PhiH2'*yCov2^(-1)*PhiH2)^(-1)));
err = [p.i err1 err2];
xCov = [diag(p.P_cov) diag(P_prior) diag(P_post)];

% 
% % Compute posterior xhat covariance
% HH   = [PhiH(:,1:10), PhiH(:,end)];
% Hbar = HH;
% % zerocol = ~any(HH,1); % assigns 0 if col of 0s is present
% % Hbar(:,zerocol) = [];   % removes cols of 0s
% % zerorow = ~any(Hbar,2); % assigns 0 if row of 0s is present
% % Hbar(zerorow,:) = [];   % removes rows of 0s
% Rbar = eye(size(Hbar,1))./p.sig_y^2;
% P_post = (Hbar'*Rbar*Hbar + P_prior^(-1))^(-1);
% 
% % Compute GDOP
% H_os2 = [H, H_offset];
% R2 = p.sig_y^2.*eye(num); % noise covariance
% warning off
% GDOP = sqrt(trace((H_os2'*R2^(-1)*H_os2)^(-1)));
% warning on
end