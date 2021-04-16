function GDOP = computeGDOP(p,H,b,num)

Phi  = diag(b(1:num));
PhiH = Phi*H(1:num,1:3);
yCov = p.sig_y^2.*eye(num); % noise covariance
GDOP = sqrt(trace((PhiH'*yCov^(-1)*PhiH)^(-1)));

end

