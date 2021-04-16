function [output] = lts(H,y,g)
% This function compute Least Trimmed Squares regression.
% Author: Ashim Neupane, email: aneup001@ucr.edu
% Date: 8-31-2020
%
% INPUT :   y - (mx1) response variable  
%           H - (mxn) predictor variable
%           g - LTS coverage
% OUTPUT:  slope - x estimate
%          flag  - selection vector

[m,n] = size(H);
if nargin == 2
%     g = floor(m+n+1)/2) + 1;
    g = floor(0.75*m) + 1;
elseif nargin == 3
    if g <= floor(m/2)
        warning('LTS: g must be more than half the number of measurements.')
        return
    elseif g > m
        warning('LTS: g must be less or equal to total number of measurements.')
        return
    end
end
% OLS  = @(y,H) (H'*H)^(-1)*H'*y; % Ordinary Least Squres
res  = @(y,H,xhat) y - H*xhat;  % compute residual
xold = OrdinaryLS(y,H); % compute initial estimate
Qold = sum(res(y,H,xold).^2);   % Cost function
dQ   = inf;  % Cost difference between successive iterations
Tol  = 1e-4; % tolerance
iter = 0;
while dQ >= Tol
    residual = res(y,H,xold);           % y - yhat
    ressquar = residual.^2;             % squared residual
    [~,idx]  = sort(ressquar,'ascend'); % index in ascending order
    short_idx = idx(1:g);               % index truncated to g entries
    short_y = y(short_idx);             
    short_H = H(short_idx,:);           
    xnew = OrdinaryLS(short_y, short_H); % OLS on selected data
    xlen = size(xnew,1);
    if (n-xlen) > 0 
        % pad xnew with 0 for every col of H removed
        xnew = padarray(xnew,n-xlen,0,'post');
    end
    Qnew = sum(res(short_y,short_H,xnew).^2);
    dQ   = abs(Qold - Qnew);
    xold = xnew;
    Qold = Qnew;    
    iter = iter + 1;
end
flag = zeros(m,1);
flag(short_idx) = 1;
output.flag  = flag;        % selection vector, 1 = selected, 0 = not
output.slope = xnew;        % x estimate
output.total_iter = iter;   % total iterations at convergence
end

function slope = OrdinaryLS(y,H)
% Compute ordinary least squares.

% check if H has any column of 0
zerocol = ~any(H,1); % check for col of 0
H(:,zerocol) = [];   % removes col of 0
% check if H has any row of 0
zerorow = ~any(H,2); % check for col of 0
H(zerorow,:) = [];   % removes col of 0
if sum(zerorow) ~= 0
    y(zerorow) = []; % removes row of 0
end

% QR is mathematically stable 
% compared to (H'*H)^(-1)H'*y 
[Q,R] = qr(H);
slope = R\(Q'*y);
end