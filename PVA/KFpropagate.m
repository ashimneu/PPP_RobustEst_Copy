function [xbatch,Pbatch] = KFpropagate(p,x_old,P_old)
% Propagate a PVA model of a moving rover

numPropSteps = p.numPropSteps;
xbatch = zeros(p.ns,numPropSteps);
Pbatch = zeros(p.ns,p.ns,numPropSteps);

for tau = 1:1:numPropSteps
    x_new = p.F*x_old;
    P_new = p.F*P_old*p.F' + p.Gam*p.Q*p.Gam';
    
    xbatch(:,tau)   = x_new;
    Pbatch(:,:,tau) = P_new;
    
    x_old = x_new;
    P_old = P_new;
end
    
end
%--------------------------------------------------------------------------