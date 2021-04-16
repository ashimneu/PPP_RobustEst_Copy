function [pos,clock_bias,res,GDOP,nsv,dnsv,nprior,dnprior,postxhat,postxhatCov] = ...
    RAPSlinear(p,xk,Pcov,H_offset,s_pos_ecef,y,eps)
% xk is the prior state
m = length(y); % total number of measurements
H = zeros(m,4);
R = zeros(m,1);
r = zeros(m,1);
off = zeros(m,1);
% Pcov = diag(Pcov);

%-----------------------%      
for j=1:m
    R(j)=norm(s_pos_ecef(:,j)-xk(1:3));
    V= (xk(1:3)-s_pos_ecef(:,j))'/R(j);
    H(j,:)=[V 1];  
    r(j) = R(j)+sagnac(p,s_pos_ecef(:,j),xk);
    if ~isempty(H_offset)
        ind = find(H_offset(j,:)==1);
        if ~isempty(ind)
            off(j) = xk(4+ind);
        end
    end
end
res = y- r - xk(4)-off;
H_os = [H,H_offset];
%-----------------------%
yCov  = p.sig_y^2.*eye(m);
E_R   = chol(yCov^(-1));
E_P   = chol(Pcov^(-1));
n     = size(H_os,2);
mu_x  = zeros(n,1);

eps_pos_x = eps; eps_pos_y = eps; eps_pos_z = eps;
Pu_pos   = [eps_pos_x^2 eps_pos_y^2 eps_pos_z^2];
eps_vel_x = geteps(0.5774,0.7); eps_vel_y = eps_vel_x; eps_vel_z = eps_vel_x;
Pu_vel   = [eps_vel_x^2 eps_vel_y^2 eps_vel_z^2];
eps_acl_x = geteps(5,0.7); eps_acl_y = eps_acl_x; eps_acl_z = eps_acl_x;
Pu_acl   = [eps_acl_x^2 eps_acl_y^2 eps_acl_z^2];
Pu_pva   = [Pu_pos, Pu_vel, Pu_acl];
Pu_clk   = geteps(3.5e+05,0.7);
Pu_drift = geteps(10,0.7);

Pu = diag([Pu_pva, Pu_clk, Pu_drift]);
J_lwr_bound = Pu^-1;

delta_x_old = mu_x;
by_old = ones(m,1);
isidentical = false;
iter = 1;
by_old_list = [];

while ~isidentical
    [delta_x1,by,~,xflag] = RAPS(res,H_os,Pcov,yCov,E_P,E_R,mu_x,J_lwr_bound,delta_x_old);
    isidentical = isequal(by_old,by);
    delta_x_old = delta_x1;
    by_old = by;
    by_old_list = insert_vector(by_old_list,by,-1);
    iter = iter + 1;
    if iter > 20
        break
    end
end
dnprior = xflag;
% dJl = diag(J_lwr_bound);
% dJp = diag(Jp);

Pbx  = eye(n);
Pby  = diag(by);
PhiH = Pby*H_os;
delta_x2 = ((PhiH'*yCov^(-1)*PhiH + Pbx'*Pcov^(-1)*Pbx)^(-1))*...
    (PhiH'*yCov^(-1)*Pby*res + Pbx'*Pcov^(-1)*Pbx*mu_x);
%-----------------------%
x_hat = xk + delta_x2;
pos = x_hat(1:3);
clock_bias = x_hat(4);
postxhat = x_hat(1:4); % Posterior xhat
% err = norm(p.truepos - pos);

% removes zero columns & rows to prevent singular matrix inversion
% when any sat constellation isn't present after measurement selection
% threshold = 0.0001;
Hbar    = PhiH;
zerocol = ~any(PhiH,1); % assigns 0 if col of 0s is present
Hbar(:,zerocol) = [];   % removes cols of 0s
zerorow = ~any(Hbar,2); % assigns 0 if row of 0s is present
Hbar(zerorow,:) = [];   % removes rows of 0s

% Compute posterior xhat covariance
yCovinvsub = eye(size(Hbar,1))./p.sig_y^2;
Pcovsub = Pcov(1:3,1:3);
if isempty(Hbar)
    postxhatCov = []; % Posterior xhat Covariance 
else
    Hsub = Hbar(:,1:3);
    postxhatCov = (Hsub'*yCovinvsub*Hsub + Pcovsub^(-1))^(-1); % Posterior xhat Covariance
end

% compute GDOP
yCovbar = p.sig_y^2.*eye(size(Hbar,1));
warning off
GDOP    = sqrt(trace((Hbar'*yCovbar^(-1)*Hbar)^(-1)));
warning on

nsv     = sum(by); % number of measurements used
dnsv    = m - nsv; % number of measurements discarded
nprior  = n; % number of priors used
% dnprior = 0; % number of priors discarded % uncomment this. temporarily
% commented.
end
