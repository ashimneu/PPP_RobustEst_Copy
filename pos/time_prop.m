function [p] = time_prop(p,dt)

lamd_a = 0.1;
lamd_c = 1.0;
a1 = exp(-lamd_a*dt);
a2 = (1-a1)/lamd_a;
a3 = (lamd_a*dt-1+a1)/(lamd_a^2);

b1 = 1;
b2 = dt;

phi_v = [eye(3,3),dt*eye(3,3),a3*eye(3,3);
         zeros(3,3),eye(3,3),a2*eye(3,3);
         zeros(3,3),zeros(3,3),a1*eye(3,3)];
Gam_v = [sqrt((dt^5)/20)*eye(3,3), zeros(3,3), zeros(3,3);
         zeros(3,3), sqrt((dt^3)/3)*eye(3,3), zeros(3,3);
         zeros(3,3), zeros(3,3), sqrt(dt)*eye(3,3)];
Q_v = [zeros(6,6), zeros(6,3);
       zeros(3,6), p.pva_siga^2*eye(3,3)];
   
phi_c = [eye(p.sys_num),[b2;0;0];
         zeros(1,p.sys_num),b1];
Gam_c = zeros(p.sys_num+1,p.sys_num+1);
Gam_c(1,1) = sqrt((dt^3)/3);
Gam_c(end,end) = sqrt(dt);
% Gam_c = [sqrt((dt^3)/3)*eye(p.sys_num),zeros(3,1);
%          zeros(1,p.sys_num),sqrt(dt)];
% Q_c = (p.pva_sigc^2)*eye(p.sys_num+1);
Q_c = zeros(p.sys_num+1,p.sys_num+1);
Q_c(1,1) = p.pva_sigc^2;
Q_c(end,end) = p.pva_sigdc^2;
Q_c(2,2) = p.pva_ISB_E^2; 
Q_c(3,3) = p.pva_ISB_B^2;
phi = [phi_v, zeros(9,p.sys_num+1);
       zeros(p.sys_num+1,9),phi_c];
Gam = [Gam_v, zeros(9,p.sys_num+1);
       zeros(p.sys_num+1,9),Gam_c];
Q = [Q_v, zeros(9,p.sys_num+1);
     zeros(p.sys_num+1,9),Q_c];
x = p.state_PVA;
x_new = phi*x;
p.state_PVA = x_new;
P_new = phi*p.P_cov*phi'+Gam*Q*Gam';
p.P_cov = P_new;

