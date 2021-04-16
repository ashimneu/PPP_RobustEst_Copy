function output = func_Q_v(sigma_pos,sigma_vel,sigma_acl,n)
    % OUTPUT - rover process noise covariance matrix associated to position,
    %          velocity, & acceleration.
    % INPUT  - sigma_acl - process noise for rover acceleration
    %                    n - dimension of rover's position/vel/accl state vector
    One       = ones(n,1);
%     sigma_pos = 0; % process noise for rover position
%     sigma_vel = 0; % process noise for rover velocity
    output    = diag([sigma_pos^2.*One; sigma_vel^2.*One;...
                sigma_acl^2.*One]);
end