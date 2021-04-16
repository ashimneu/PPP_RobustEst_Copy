function [p,res] = userpos(p,cpt)
% This is solver for computing user position, receiver clock
% bias, satellite system offset.
% Measurement selection applied
% Input: 
%       s_pos_ecef: 3-by-N Satellite position in ECEF frame.
%       x0 : 3-by-1 initial interative coordinates in ECEF frame.
%       y: m-by-1 Corrected pseudorange.
%
% Output:
%       
%       
%       
%       

%-------------------%
% Initialize

ind = find(cpt.num_sv([1,3,4]) ~= 0);
xk = [p.state_PVA(1:6);p.state_PVA(9+ind);p.state_PVA(13)];
xp = [p.state_PVA(1:3);p.state_PVA(9+ind)];
xv = [p.state_PVA(4:6);p.state_PVA(13)];
% [H_offset,x_offset] = sys_offset(cpt.num_sv);
%------------------%
% [pos,clock_bias,res] = LSsolver(p,xk,H_offset,cpt);
[xk,res] = LSsolver_PVA(p,xk,cpt);
% [xk,res] = LSsolver_test(p,xp,xv,cpt);
p.state_PVA([1:6,9+ind,13]) = xk;


% switch p.select
%     case 0
%         [pos,clock_bias,res] = LSsolver(p,xk,H_offset,s_pos_ecef,y);
%     case 1
%         [pos,clock_bias,res,cost] = TDsolver(p,xk,H_offset,s_pos_ecef,y);
%     case 2
%         [pos,clock_bias,res,cost] = LSSsolver(p,xk,H_offset,s_pos_ecef,y);
%     case 3
%         [pos,clock_bias,res,cost] = MSsolver(p,xk,H_offset,s_pos_ecef,y);
%     case 4
%         [pos,clock_bias,res,cost] = LTSsolver(p,xk,H_offset,s_pos_ecef,y);
% end

end
