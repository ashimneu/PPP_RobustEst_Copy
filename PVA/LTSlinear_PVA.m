function [pos,clock_bias,res,GDOP,nsv,dnsv,nprior,dnprior,x_post,P_post] = ...
    LTSlinear_PVA(p,xk,Pcov,H_offset,s_pos_ecef,y,x_prior,P_prior_in)
% xk is the prior state
m = length(y); % total number of measurements
H = zeros(m,4);
R = zeros(m,1);
r = zeros(m,1);
off = zeros(m,1);
xk(1:3) = x_prior(1:3);
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
%--------------------------------------------------------------------------
H_pos = H(:,1:3);                       % H for position component of state vector
H_vel = zeros(m,3);                     % H for velocity component
H_acl = zeros(m,3);                     % H for acceleration component
H_c   = [H(:,4), H_offset, zeros(m,1)]; % H for clock bias & drift components
H_os  = [H_pos, H_vel, H_acl, H_c];
[x_prior,P_prior] = adjust_entries(p,x_prior,P_prior_in,1);
%-----------------------%
R    = p.sig_y^2.*eye(m);
E_R  = chol(R^(-1));
E_P  = chol(P_prior^(-1));
n    = size(H_os,2);
mu_x = zeros(n,1);
A    = [E_R*H_os; E_P];
c    = [E_R*res; E_P*mu_x];

Abar    = A;
zerocol = ~any(A,1); % assigns 0 if col of 0s is present
Abar(:,zerocol) = [];  % removes cols of 0s

Option = 1; 
switch Option
    case 1
        [rew,~] = ltsregres(Abar,c,'plots',0,'intercept',0);
    case 2
        h = floor((m+n+n+1)/2)+1; % for highest breakdown point
        rew = lts(A,c,h);
end 

b_LTS   = rew.flag; % binary vector corresponding to measurement selection
by      = b_LTS(1:m);
bx      = b_LTS(m+1:end);
nsv     = sum(by); % number of measurements used
dnsv    = m - nsv; % number of measurements discarded
nprior  = sum(bx); % number of priors used
dnprior = n - nprior; % number of priors discarded
% mslc_count = [nsv; dnsv; nprior; dnprior]; % msr selection count

Pby = diag(by);
Pbx = diag(bx);
PhiH = Pby*H_os;
delta_x = ((PhiH'*R^(-1)*PhiH+ Pbx'*P_prior^(-1)*Pbx)^(-1))*...
   (PhiH'*R^(-1)*Pby*res + Pbx'*P_prior^(-1)*Pbx*mu_x);
%-----------------------%
x_hat = x_prior + delta_x;
pos   = x_hat(1:3);
clock_bias = x_hat(4);
x_post = x_hat; % Posterior xhat
% err = norm(p.grdpos - pos); % error check for debugging

[x_post,P_prior,Pbx] = adjust_entries(p,x_post,P_prior,2,Pbx);

% % removes zero columns & rows to prevent singular matrix inversion
% % when any sat constellation isn't present after measurement selection
% Hbar    = PhiH;
% zerocol = ~any(PhiH,1); % assigns 0 if col of 0s is present
% Hbar(:,zerocol) = [];   % removes cols of 0s
% zerorow = ~any(Hbar,2); % assigns 0 if row of 0s is present
% Hbar(zerorow,:) = [];   % removes rows of 0s

% Compute posterior xhat covariance
Hbar    = [PhiH(:,1:10), PhiH(:,end)];
Rinvbar = eye(size(Hbar,1))./p.sig_y^2;
J_post  = Hbar'*Rinvbar*Hbar + Pbx'*P_prior^(-1)*Pbx;
% J_post2  = remove0colrow(J_post,'column');
P_post 	= J_post2^(-1);

% % Compute posterior xhat covariance
% yCovinvsub = eye(size(Hbar,1))./p.sig_y^2;
% Pcovsub = Pcov(1:3,1:3);
% Hsub = Hbar(:,1:3);
% postxhatCov = (Hsub'*yCovinvsub*Hsub + Pcovsub^(-1))^(-1); % Posterior Covariance

% Compute GDOP
H_os2 = [H, H_offset];
R2 = p.sig_y^2.*eye(m); % noise covariance
GDOP = sqrt(trace((H_os2'*R2^(-1)*H_os2)^(-1)));

end