function delta = sagnac_v(p,sat,x)
% Compute Sagnac correction/earth rotate correction
% sat,x: {x,y,z,vx,vy,vz} in ECEF
% Reference: https://escholarship.org/uc/item/1bf6w7j5
% (omge/c)*(vs_y * p_x + ps_y*v_x - vs_x*p_y-ps_x*v_y)
delta = p.omge*(sat(5)*x(1)+sat(2)*x(4)-sat(4)*x(2)-sat(1)*x(5))/p.c;
end