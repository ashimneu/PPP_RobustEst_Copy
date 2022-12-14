clear all
format compact

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters of the simulation:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tau = 0.1;  % simulation time step
T=1;        % sampling interval for sensors
N_sm_stps = ceil(T/tau);
tau = T/N_sm_stps;  % ensure integer number of sim steps in each T
T_end = 240; %  duration of simulation

% covariance matrix of position measurements
R_Pos = 1;
sigma_P = sqrt(R_Pos);
m  = 25; % number of measurements

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%real world simulation:
% This creates:
% X_gt - the actual (ground truth) vehicle trajectory. This is only for analysis of
% performance and should not be used in the state esimator. This is what you
% are trying to estimate.
% t_gt - times of the ground truth state
%
% z_pos - Measurements linearly related to position at t = k T.
% t_z   - times of z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
[X_gt, z_pos, t_gt, t_z] = rws(tau, T, T_end, sigma_P, m);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ToDo 2: Add outliers later, after state estimator is properly working. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ToDo 1: Implement and run the state estimator(s) here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%% x-y plot of true trajectory
figure(1)
clf
plot(X_gt(1,:),X_gt(2,:));
hold on
grid on
% add axis labels
xlabel ('x (m)')
ylabel ('y (m)')


figure(2)
clf
ylabels =[' North, N, m '
          ' East, E, m  '
          ' Down, D, m  '
          ' Vel. N, mps '
          ' Vel. E, mps '
          ' Vel. D, mps '
          'Acc. N, mpsps'
          'Acc. E, mpsps'
          'Acc. D, mpsps'];
for i=1:9
    ax(i) = subplot(3,3,i);
    plot(t_gt,X_gt(i,:))
    grid on
    xlabel('Time, t, s')
    ylabel(ylabels(i,:))
end
linkaxes(ax,'x')

 
figure(3)
clf
plot(t_z,z_pos,'.')
grid on
xlabel('Time, t, sec.')
ylabel('measurements, m')
 