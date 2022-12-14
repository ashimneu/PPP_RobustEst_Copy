function [X, y,  t_gt, t_msr] = rws(dt, T, Tf, sigma_P, m)
% Implements a point particle PVA linear simulation. 
% Computes m measurements linearly related to the position
% each with GRV meassurement noise with std sigma_P
t_gt = 0:dt:Tf; % simulation times for X_gt
t_msr= 0:T:Tf;  % times of measurements y
H    = [rand(m,3),zeros(m,6)];   % measurement matrix
for i=1:m
    H(i,1:3) = H(i,1:3)/norm(H(i,1:3));  % make it a unit vector
end
RE = 0.7;
IM = 0.5;
p1 = RE + IM*sqrt(-1);
p2 = RE - IM*sqrt(-1);
K(1) = p1 * p2;
K(2) = (p1+p2)-1;

% Preallocate the matrices (this speeds up the code, by allocating memory 
% just once, and not every time a new column is added to the matrix. 
% It can make a huge difference in runtime when the number of 
% timesteps is large)
Nx = length(t_gt);
X  = zeros(9,Nx);
u  = zeros(3,Nx);
Ny = length(t_msr);
y  = zeros(m,Ny);
 
I3 = eye(3,3);
Z3 = zeros(3,3);
Phi = [I3 I3*dt 0.5*I3*dt*dt 
       Z3  I3     I3*dt
       Z3  Z3     I3];
G   = [Z3; Z3; I3];
omga = 2*pi/Tf;
z_ampl = 2;
A = z_ampl*omga; 
nmsr = 0;
% first measurement at time 0
nmsr = nmsr + 1;
y(:,nmsr) = H * X(:,i) + sigma_P*randn(m,1); % noise corrupted measurement
t_msr(nmsr) = t_gt(1);
for i=2:Nx
    % define commanded speed to trace out a box
    Vc = zeros(3,1);
    Vc(3) = A * sin(omga * t_gt(i));
    delt = rem(t_gt(i),120)
    if delt < 10
        Vc(1) = delt;
    elseif delt < 20
        Vc(1) = 10;
    elseif delt < 30
        Vc(1) = 30 - delt;
    elseif delt < 40
        Vc(2) = delt-30;
    elseif delt < 50
        Vc(2) = 10;
    elseif delt < 60
        Vc(2) = 60 - delt;
    elseif delt < 70
        Vc(1) = 60 - delt;
    elseif delt < 80
        Vc(1) = -10;
    elseif delt < 90
        Vc(1) = delt - 90;
    elseif delt < 100
        Vc(2) = 90 - delt;
    elseif delt < 110
        Vc(2) = -10;
    elseif delt < 120
        Vc(2) = delt - 120;
    end
    u = K(1) * (Vc - X(4:6,i-1)) - K(2) * X(7:9,i-1); % speed control law
    X(:,i) = Phi * X(:,i-1) + G * u;                  % dynamic model: PVA
    if abs(t_gt(i) - t_msr(nmsr+1)) < dt/2            % find measurement times
        nmsr = nmsr + 1;
        y(:,nmsr) = H * X(:,i) + sigma_P*randn(m,1);  % noise corrupted measurement
        t_msr(nmsr) = t_gt(i);
    end
end
nmsr




