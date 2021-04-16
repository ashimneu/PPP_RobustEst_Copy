function [pos,msr_res,GDOP,nsv,dnsv,nprior,dnprior,x_post,P_post,dx,Jydiag,by,H_pos_vel] = ...
    MSHuberlinear_PVA_useallprior(p,cpt,grdpos,x_prior,P_prior,HuberConst)
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
yCov = blkdiag(p.sig_y^2.*eye(num),p.sig_y_dop^2.*eye(num)); % noise covariance
E_R  = chol(yCov^(-1));
E_P  = chol(Pcov^(-1));
n     = numel(xk);
mu_x  = zeros(n,1);
A     = [E_R*H_os; E_P];
c     = [E_R*[res_R;res_D]; E_P*mu_x];

flag = [NaN(2*num,1); ones(n,1)];
[dx,Wsquared] = Mestimator(c,A,HuberConst,flag);
w = sqrt(diag(Wsquared)); % vector of weights
%-----------------------%
xk  = xk + dx;
x_post = xk;
pos = xk(1:3);
err1 = norm(grdpos - pos);
end

%---------------------------
wy = w(1:2*num); by = wy;
wx = w(2*num+1:end);
nonzero_wy = find(wy); % find non-zero weights
nonzero_wx = find(wx); % find non-zero weights
nsv  = numel(nonzero_wy); % number of measurements used
dnsv = num - nsv; % number of measurements discarded
nprior  = size(nonzero_wx,1);  % number of priors used
dnprior = size(wx,1) - nprior; % number of priors discarded

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

% for GDOP computation, truncate A & weights to remove portion
% corresponding to prior & keep that corresponding to measurement.
A_trunc = A(1:2*num,:); % truncated A matrix
Uy = diag(wy);          % weight matrix truncated
aug_A = Uy * A_trunc;   % 

% Remark: A_aug is likely to have 0's column after measurement selection by
% weights_trunc when all satellites observation of certain constellation
% are given 0 weight by Mestimator. This will lead to singular matrix 
% inversion in GDOP calculation. To prevent this, 0's column need to be 
% checked for and removed.

% remove 0's columns from A_aug (if present)
aug_Abar = aug_A;
zerocol  = ~any(aug_A,1); % index of 0's column
aug_Abar(:,zerocol) = []; % removes 0's column

% compute GDOP
GDOP = sqrt(trace(inv(p.sig_y^-2*(aug_Abar'*aug_Abar))));

%--------------------------%
% Compute residual
E_R = sqrtm(yCov^(-1));
E_P = sqrtm(P_prior^(-1));
Ab  = [E_R*Uy*H_os; E_P];
cb  = [E_R*Uy*[res_R;res_D]; E_P*mu_x];
residual = cb - Ab*dx ;
msr_res  = residual(1:2*num);

% Compute posterior xhat covariance
Hbar   = Uy*H_os;
H_pos_vel = Hbar(:,1:6);
Rbar_R = eye(num)./p.sig_y^2; 
Rbar_D = eye(num)./p.sig_y_dop^2;
Rbar   = blkdiag(Rbar_R,Rbar_D);
P_post = (Hbar'*Rbar*Hbar + P_prior^(-1))^(-1);

% Compute posterior information (obtained from observations only)
Jydiag = diag(Hbar'*Rbar*Hbar);

end
%%%%=======================================================================
function [xEstimate,Weights] = Mestimator(y,H,const,flag)
    m = size(H,1);
    if nargin == 3
        flag = NaN(m,1);
    end
    Tol = 1e-4;    % IRWLS iteration tolerance
    Totaliter = 0; % total iterations before convergence
    K = 1.345;     % tuning constant for Huber objective function
    
    if nargin == 3
        K = const;
    end    
    
    % Caclulate Ordinary LS
    OLS  = @(y,H) (H'*H)\H' * y; % Ordinary Least Squres
    xEst_old = OLS(y,H); % compute inital estimate
    diff = inf;   
    iterCount = 0;
    
    %Iteratively Reweighted Least Squares method 
    while diff >= Tol && iterCount < 1000       
        res    = y - H*xEst_old;
        sigEst = calcSigEst(res);
        u         = res./sigEst;        
        [xEst,W]  = calcWLS(u,y,H,K,flag);
        diff      = norm(xEst - xEst_old);
        xEst_old  = xEst;
        iterCount = iterCount + 1; 
    end  
    xEstimate = xEst;
    Weights   = W;
    Totaliter = iterCount;
end
%--------------------------------------------------------------------------
function [xEst,W] = calcWLS(u,y,H,const,flag)
    % const - normalizing constant for Huber objective function
    % calculates Weighted Least Squares
    m = size(H,1);
    n = size(H,2);
    w = zeros(m,1);
    for i = 1:1:m
        flg = flag(i);
        if isnan(flg)
            w(i) = calcWeight(u(i),const);  % get weights for u_i
        else
            w(i) = flg;
        end
    end
    W    = diag(w);    
%     xEst = (H'*W*H)^(-1)*(H'*W*y);        % Weighted LS
    V = sqrt(W);
    Hbar = V*H;
    zerocol = ~any(Hbar,1); % assigns 0 if col of 0s is present
    Hbar(:,zerocol) = [];   % removes cols of 0s    
    xEst = inv(Hbar'*Hbar)*(Hbar'*V*y);
    xlen = size(xEst,1);
    if (n-xlen) > 0 
        % pad xEst with a 0 for every satellite constellation 
        % measurements ignored
        xEst = padarray(xEst,n-xlen,0,'post');
    end
end
%--------------------------------------------------------------------------
function SigEst = calcSigEst(res)
% calculates sigma_hat, estimate of scale for Regression M-Estimator

    MAD     = median(abs(res-median(res))); % Max Absolute Deviation (MAD)
    factor  = 0.675;      % MAD normalization factor
    SigEst  = MAD/factor; % Normalized MAD    
end
%--------------------------------------------------------------------------
function weight = calcWeight(u,const)
    % Caluclates weights for Weighted Least Squares
    % const - normalizing constant
    % u     - res/sig 

    if abs(u) <= const
        weight = 1;
    else
        weight = const/abs(u);
    end
end
%--------------------------------------------------------------------------