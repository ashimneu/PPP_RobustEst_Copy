function Epsilon = geteps(alpha_beta)
% Computes performance specification parameter, epsilon, for RAPS. 
% See eqn. (28) of [1].
% [1] - E. Aghapour, F. Rahman, J. A. Farrell, Outlier Accomodation in
% Nonlinear State Estimation: A Risk-Averse Performance-Specified Approach
% INPUT:    alpha - [meters] Real number
%           beta -  error percentage, value should be between 0 & 100%,
    alpha = alpha_beta(:,1);
    beta  = alpha_beta(:,2)./100;
    Epsilon = alpha./(sqrt(2)*erfinv(beta));
end

