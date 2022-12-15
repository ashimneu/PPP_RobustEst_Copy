function [x_post,by,augcost,exitflag] = RAPS(y,H,P,R,J_l,x_prior)
% Finds all binary vectors that statisfy performance constraint of RAPS. 
% Next, uses LS approach to estimate state x which corresponds to the
% lowest cost function.
% OUTPUT:   x_post   - posterior state vector estimate
%           by       - measurement selection vector (binary)
%           augcost  - augmented cost for RAPS B&B
%           exitflag - see MATLAB function: intlinprog for description
% INPUT:    y    - measurements
%           H    - measurement matrix
%           Pcov - Prior Covariance matrix
%           yCov - Measurement Covariance matrix
%           J_l  - Information Matrix Lower Bound
%           x_prior - prior state vector estimate
    
    Jpminus = P^-1; % state prior info. matrix
    Jrminus = R^-1; % measurement info. matrix
    E_P = chol(Jpminus); % cholesky decomp. of state prior info. matrix
    E_R = chol(Jrminus); % cholesky decomp. of measurement info. matrix

    [m,n] = size(H);
    lowerbound = zeros(m,1); % lower bound on measurement selection vector
    upperbound = ones(m,1);  % upper bound on measurement selection vector
    
    diagR  = diag(R);               % diagonal entries of measurement covariance
    ieqLHS = -((H').^2)./repmat(diagR',n,1); % inequality constraint LHS matrix
    ieqRHS = diag(Jpminus - J_l);   % inequality constraint RHS vector   
    cost   = (y-H*x_prior).^2 .* diag(Jrminus) ; % cost function
    intcon = 1:1:m;                 % entries of by vector that take only integer value
    option = optimoptions(@intlinprog,'display','off'); % for output supression
    
    % solve for optimal selection vector (uses Branch & Bound search)
    [by,~,exitflag,~] = intlinprog(cost,intcon,ieqLHS,ieqRHS,[],[],lowerbound,upperbound,option);
    

    bx = ones(n,1); % prior state selection vector (binary) 
    if exitflag == 1
        % feasible solution is found
        [x_post, augcost]  = MAP(by,y,H,E_R,E_P,x_prior);
        J_out   = calcJb(by,H,R,Jpminus);        
    else
        % feasible solution not found
        fprintf('\n Exitflag = %1.0f. Feasible solution not found.\n',exitflag)
        fprintf('Check exitflag in function "intlinprog" for more details.\n')
        
        fprintf('Selecting all measurements.\n')
        by = ones(m,1);
        x_post  = MAP(by,y,H,E_R,E_P,x_prior);
        augcost = Inf;
        J_out   = calcJb(by,H,R,Jpminus); %#ok<*NASGU> 
    end
        
end
%--------------------------------------------------------------------------
function [x_post,aug_cost] = MAP(by,y,H,E_R,E_P,x_prior)
    % Maximum A Posteriori state estimate and the augmented cost function
    % when doing measurement selection
    % INPUT: by  - measurement selection vector
    %        E_P - sqrt(P^-1)
    %        E_R - sqrt(R^-1)

    if isempty(by)
        error('Input arg "by" is empty.');
    else
        Phiby = diag(by);
        A = [E_R * Phiby * H; E_P];
        c = [E_R * Phiby * y; E_P * x_prior];
        x_post = (A'*A)^-1*A'*c; % posterior state estimate
        aug_cost = norm(A*x_post-c)^2;
    end

end
%--------------------------------------------------------------------------
function J = calcJb(by,H,R,Jpminus)
    % Calculates posterior state information matrix when doing measurement selection
    % INPUT: by - measurement selection vector
    %        Jpminus - State information matrix prior, Jpminus = Pminus^-1
    
    Pby  = diag(by); % Phi(by)                        
    PhiH = Pby * H;
    J    = PhiH' * R^-1 * PhiH + Jpminus;
end