function [xhat,by,Augcost,exitflag] = RAPS(y,H,Pcov,yCov,E_P,E_R,mu_x,Jl,xhat_old,iter,p)
% Finds all binary vectors that statisfy performance constraint from RAPS. 
% Next, uses LS approach to find x corresponding to lowest cost and 
% associated binary vector.
    if nargin == 8
        xhat_old = mu_x;
    end
    [m,n] = size(H);
    dR = diag(yCov); % diagonal of Covariance matrix of y
    lowerbound = zeros(m,1); % binary vector lower bound           
    upperbound = ones(m,1);  % binary vector upper bound 
    Jminus  = Pcov^-1;
    ieqLHS  = -((H').^2)./repmat(dR',n,1); % inequality constraint LHS matrix
    ieqRHS  = diag(Jminus - Jl);  % inequality constraint RHS vector   
    objfunc = (y-H*xhat_old).^2 .* diag(yCov^-1) ; % objective function
    intcon  = 1:1:m; % entries of by that take only integer value
    options = optimoptions(@intlinprog,'display','off'); % for B&B output supression
    
    % solve for optimal selection vector by using (Branch & Bound search)
    [by,~,exitFlag,~] = intlinprog(objfunc,intcon,ieqLHS,ieqRHS,[],[],lowerbound,upperbound,options);
    bx = ones(n,1);

    if exitFlag == 1
        % feasible solution is found.
        xhat = MAP(bx,by,y,H,E_R,E_P,mu_x);
        Augcost = getAugcost(by,xhat,y,H,E_R,E_P,mu_x);
%         J_out = calcJb(by,H,yCov,Pcov);        
    else
        % Feasible solution not found
%         fprintf('\n exitflag = %1.0f. Feasible solution not found.\n',exitFlag)
%         fprintf(' Check "intlinprog" for descriptions on exitflag.\n')        
        by   = ones(m,1);
        xhat = MAP(bx,by,y,H,E_R,E_P,mu_x);
        Augcost = Inf;
%         J_out  = calcJb(by,H,yCov,Pcov);
    end
    
%     if sum(by(1:floor(m/2))) <= floor(m/4)
%         bysum = sum(by(1:floor(m/2)));
%     end
    
    if nargout == 4
        exitflag = exitFlag;
    end  
    
%     debug = 0;
%     if (p.i==30 && iter >1), debug = 1; end
%     if debug
%         idx = [1:6, 10:size(ieqRHS,1)]';
%         tempLHS = sum(ieqLHS,2);
%         tempLHS = tempLHS(idx,:);
%         RHS = ieqRHS(idx);
%     end
end
%--------------------------------------------------------------------------
function xhat = MAP(bx,by,y,H,E_R,E_P,mu_x)
% get standard Maximum A Posteriori estimation for 
% measurement selection using vector b.
    Pby  = diag(by);
    Pbx  = diag(bx);
    A    = [E_R*Pby*H; E_P*Pbx];
    c    = [E_R*Pby*y; E_P*Pbx*mu_x];
    xhat = (A'*A)^-1*A'*c; % MAP solution
end
%--------------------------------------------------------------------------
function aug_cost = getAugcost(by,xhat,y,H,E_R,E_P,mu_x)
    % compute augmented cost for RAPS B&B
    % by = measurement selection binary vector
    
    if isempty(by)
        warning('Input arg "by" is empty.')
        aug_cost = inf;
    else
        Phiby = diag(by);
        A = [E_R*Phiby*H; E_P];
        c = [E_R*Phiby*y; E_P*mu_x];
        aug_cost = norm(A*xhat-c)^2;
    end
end
%--------------------------------------------------------------------------
function Jb = calcJb(by,H,yCov,Pcov)
% calculates posterior information matrix for the selection vector b

    Pby  = diag(by); % Phi(by)                        
    PhiH = Pby * H;
    Jb   = PhiH' * yCov^-1 * PhiH + Pcov^-1;
end