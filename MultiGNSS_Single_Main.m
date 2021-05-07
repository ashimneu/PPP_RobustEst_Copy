% This code implement GNSS (GPS, GLO, GAL, BDS) under single frequency
% Special for CNES SSR data
% Correction type: PPP (Precise Point Positioning)
% clear all
% close all
%--------------------------------%
addpath('data')
addpath('parser')
addpath('time_compute')
addpath('eph')
addpath('pos')
addpath('corr')
addpath('PVA')
addpath('PVA/multisim')
addpath('pos/LIBRA')
addpath('PVA/outliers')
%--------------------------------%
% Pick the Data Number
initpath = 'data/';
data_num = 6;
[eph_name,obs_name,IGS_name,data_base,code_bia,Grdpos,USTEC] = datapathload(data_num,initpath);
%--------------------------------%

% Initialize parameters
[p,eph,obs] = initialization(eph_name,obs_name,Grdpos);
% obs.GLO.S1 = obs.GLO.S2; % Uncomment when GLO in data 8
% obs.BDS.P1(18,:) = 0;
% Mode setting
p.run_mode = 0;
p.post_mode = 1; %%%% 0=Standard GNSS, 1 = PPP, 2= DGNSS
p.VRS_mode = 0;
p.IGS_enable = 1;
p.double_diff = 0;
p.elev_mark  = 15*pi/180;
p.enableGPS  = 1; % Enable GPS: 1 means enable, 0 means close
p.enableGLO  = 0; % Enable GLO: 1 means enable, 0 means close
p.enableGAL  = 1; % Enable GAL: 1 means enable, 0 means close
p.enableBDS  = 1; % Enable BDS: 1 means enable, 0 means close
p.inval = 1; % Computation time interval
p.tec_tmax = 15;
p.tec_tmin = 0;
p.L2enable = 0;
%--------------------------------%
% PVA settings
Z3 = zeros(3,1); O3 = ones(3,1);
p.sys_num = 3;
p.state_PVA = zeros(3+3+3+p.sys_num+1,1); % P;V;A;tr+ISBs;dr(1X1)
p.pva_siga = 0.005; %1;
p.pva_sigc = 15;
p.pva_sigdc = 0.1; %2;
p.pva_ISB_E = 0.005; %1; 
p.pva_ISB_B = 0.005; %1;
% p.P_cov = ones(3+3+3+p.sys_num+1,1);
% p.P_cov(10:12) = p.pva_sigc^2;
p.P_cov = [0.1.*O3; 0.01.*O3; 0.01.*O3; 10^2; 0.002^2; 0.003^2; .1];
p.P_cov = diag(p.P_cov);
p.PVA_enable = 1;
% p.P_cov(end) = p.pva_sigc^2;
% p = PVA_parameters(p);
%--------------------------------%
% For Measurement selection
p.x0 = [-2430696.646;-4704190.763;3544329.081;0]; % prior state(1:3)
p.priorposcov = [2;2;2;1]; % cov([x0;clk;x_off])
p.ISBglo = 0.45; p.ISBglo_cov = 0.257^2; % Inter system bias GPS to GLO
p.ISBgal = 0.45; p.ISBgal_cov = 0.257^2; % Inter system bias GPS to GAL
p.ISBbds = 0.92; p.ISBbds_cov = 0.257^2; % Inter system bias GPS to BDS 
p.sig_y  = 0.951; % [meters] PPP pseudorange residual measurement covar
p.sig_y_dop = 0.25; % [m/s] PPP doppler residual measurement covar
% Enable Measurement selection Algorithms
p.eb_LTS  = 1;
p.eb_RAPS = 1;
p.eb_TD   = 1;
p.eb_MShb = 1;
p.eb_MStk = 1;
p.eb_LSS  = 0;
% Declare Algo. parameters
p.LTSOption = [2 3]; % 1 =default (check LTSlinear.m)
p.RAPSEps = [0.5 99.5; 0.4 99.9]; %; 1 99.5]; % alpha & beta
p.TDLambda  = [2.5 3 3.5];
p.MShbConst = [1.345];
p.MStkConst = [4.685];
p.LSSLambda = [1];
% Initialize parameters
p.LTSn  = numel(p.LTSOption);
p.RAPSn = size(p.RAPSEps,1);
p.TDn   = numel(p.TDLambda);
p.MShbn = numel(p.MShbConst);
p.MStkn = numel(p.MStkConst);
p.LSSn = numel(p.LSSLambda);
%--------------------------------%
% Declare & Initialize PVA parameters
Z3 = zeros(3,1); O3 = ones(3,1);
p.n  = 3; % dimension of rover position/vel/acc state vector 
p.Ms = 0; % multipath states dimension
p.ns = 13 + p.Ms;   % KF state vector dimension   
p.numPropSteps = 10; % propagation steps
p.T = 1/p.numPropSteps; % [sec] time between consecutive propagation steps
p.lam_a = 0.1;
p.lam_cdrift = 1.0;
% process noise std
p.sig_p = 0; 
p.sig_v = 0; 
p.sig_a = 0.003; % 0.0002;
p.sig_cbias  = 2; % 0.5
p.sig_cdrift = 0.1;
p.sig_ISB_E = 0.0002; 
p.sig_ISB_B = 0.004; 
% initial state
p.x_v_prior = [p.Grdpos.pos(:,1); Z3; Z3];
p.x_c_prior = [350542.5; 0.14; 10; 142];
% initial covar
p.pva_cov_prior = [0.1^2.*O3; 0.01^2.*O3; 0.002^2.*O3];
p.clk_cov_prior = [10^2; 0.01^2; 0.005^2; .1];
p = initPVAparams(p);
%--------------------------------%
p.eb_outlier  = 1;
p.genOutlier  = 1;
p.saveOutlier = 0;
p.multisim_outliervar = "mean"; % mean/width/count
p = initOutlierparam(p);
%--------------------------------%
[p,obs] = load_PPP_corr(p,data_base,IGS_name,eph,obs,USTEC,code_bia);
%-------------%
outputcell = compute_gnss_ecef_multisim(p,eph,obs);
output = outputcell{1};
save('outputcell_May_outliermean_1thru20.mat','outputcell','-v7.3');
%%
clc
opt.movingrover = 0; % For rover: 0 = stationary, 1 = moving
opt.LTSn_i  = 2;
opt.RAPSn_i = 2; 
opt.TDn_i   = 3;
opt.MShbn_i = 1;
opt.MStkn_i = 1;
opt.LSSn_i  = 0;
opt.axes2link ='x'; 
opt.ax = [];
opt.frame = "ned";
% opt = ploterrgdopscat(output,"kf",opt);
% opt = plotTraj(output,"kf",opt);
% opt = plotTraj(output,"raps",opt);

% opt = plotErrorStd(output,opt);

metrics(output,opt)
% CDF_curves = ploterrcdf(output,opt);

compare_err_std(outputcell,opt,"mean","std")

linkaxes2(opt);
% compare_err_std(outputcell,opt,"count","std")

%%
if false
opt.LTSn_i = 1;
startTime = 50;
endTime = output.N;
Timeline = startTime:1:endTime; 
fprintf('\nLS error Std. = %2.3f \n',nanstd(output.err_LS(1,Timeline)))
fprintf('RAPS error Std. = %2.3f \n',nanstd(output.err_RAPS(opt.RAPSn_i,Timeline)))

end
%%
if false
msr_res_nan = output.resLS{1};
msr_res = msr_res_nan(any(msr_res_nan,2),:); % get rows that aren't all nan
epoch_count = size(msr_res,2);
res_pse = []; res_dop = [];
for j = 1:epoch_count
    jth_res_nan = msr_res(:,j);
    jth_res = jth_res_nan(any(jth_res_nan,2));
    num = size(jth_res,1)/2;
    res_pse = [res_pse; jth_res(1:num)]; 
    res_dop = [res_dop; jth_res(num+1:end)];
end
fprintf('Pseudorange residual Std. = %2.3f \n',nanstd(res_pse))
fprintf('Pseudorange residual MAD = %2.3f \n',mad(res_pse,1))
fprintf('Doppler residual Std. = %2.3f \n',nanstd(res_dop))
fprintf('Doppler residual MAD = %2.3f \n',mad(res_dop,1))
end

%%
if false
Time = output.gpst; % p.t;
figure
scatter(Time,output.err,'.')
% xtickformat('yyyy-MM-dd HH:mm:ss')
title('ECEF positioning error')
xlabel('Local time')
ylabel('Error, unit: meter');grid on
ax = [findall(gcf, 'type', 'axes')];

figure
scatter(Time,output.hor_err,'.')
% xtickformat('yyyy-MM-dd HH:mm:ss')
title('Horizontal positioning error')
xlabel('Local time')
ylabel('Error, unit: meter');grid on
ax = [ax; findall(gcf, 'type', 'axes')];

total = output.sv_num_GPS + output.sv_num_GAL + output.sv_num_BDS;

figure
scatter(Time,output.sv_num_GPS,'.')
hold on
scatter(Time,output.sv_num_GAL,'.')
hold on
scatter(Time,output.sv_num_BDS,'.')
hold on
scatter(Time,total,'.')
title('total satellites been used')
legend('GPS','GAL','BDS','Total')
xlabel('Receiver time using GPS second')
ylabel('Distance, unit: meter');grid on
ax = [ax; findall(gcf, 'type', 'axes')];
% 
figure
subplot(311)
mps2kmph = 1; % 1/(5/18);
scatter(Time,output.v_ecef(1,:)*mps2kmph,'.')
% xtickformat('yyyy-MM-dd HH:mm:ss')
% title('ECEF Velocity')
% xlabel('Local time')
% ylabel('unit: km/h');grid on
ylabel('unit: ms^{-1}');grid on

subplot(312)
scatter(Time,output.v_ecef(2,:)*mps2kmph,'.')
% xtickformat('yyyy-MM-dd HH:mm:ss')
% title('ECEF Velocity')
% xlabel('Local time')
% ylabel('unit: km/h');grid on
ylabel('unit: ms^{-1}');grid on

subplot(313)
scatter(Time,output.v_ecef(3,:)*mps2kmph,'.')
% xtickformat('yyyy-MM-dd HH:mm:ss')
% title('ECEF Velocity')
xlabel('Local time')
% ylabel('unit: km/h');grid on
ylabel('unit: ms^{-1}');grid on
sgtitle('EKF - Velocity')
ax = [ax; findall(gcf, 'type', 'axes')];
 
figure
scatter(Time,output.rover_clk,'.')
title('Local bias')
xlabel('Receiver time using GPS second');
ylabel('Clock bias, seconds');grid on
ax = [ax; findall(gcf, 'type', 'axes')];
% 
figure; hold on
for kk = 1:size(output.res_GPS,1)
scatter(Time,output.res_GPS(kk,:),'.')
end
title('Residual GPS')
xlabel('Receiver time using GPS second');
ylabel('residual, meters');grid on
ax = [ax; findall(gcf, 'type', 'axes')];

figure; hold on
for kk = 1:size(output.res_GAL,1)
scatter(Time,output.res_GAL(kk,:),'.')
end
title('Residual GAL')
xlabel('Receiver time using GPS second');
ylabel('residual, meters');grid 
ax = [ax; findall(gcf, 'type', 'axes')];

figure; hold on
for kk = 1:size(output.res_BDS,1)
scatter(Time,output.res_BDS(kk,:),'.')
end
title('Residual BDS')
xlabel('Receiver time using GPS second');
ylabel('residual, meters');grid 
ax = [ax; findall(gcf, 'type', 'axes')];

figure; hold on
for kk = 1:size(output.msr_res_GPS,1)
scatter(Time,output.msr_res_GPS(kk,:),'.')
end
title('Linear Msr Residual GPS')
xlabel('Receiver time using GPS second');
ylabel('residual, meters');grid 
ax = [ax; findall(gcf, 'type', 'axes')];


if ~isempty(ax) 
    if isfield(opt,'ax')
        opt.ax = [opt.ax; ax];
    else
        opt.ax = ax;
    end
end

linkaxes2(opt);

end

% linkaxes(opt.ax,opt.axes2link)

% figure
% scatter(p.t,output.res_GPS,'.')
% title('Residual GPS')
% xlabel('Receiver time using GPS second');
% ylabel('Clock bias, seconds');grid on

% figure
% subplot(311)
% scatter(p.t,output.ned_err(1,:),'.')
% title('North Error in NED');grid on;
% subplot(312)
% scatter(p.t,output.ned_err(2,:),'.')
% title('East Error in NED');grid on;
% subplot(313)
% scatter(p.t,output.ned_err(3,:),'.')
% title('Down Error in NED')
% xlabel('Receiver time using GPS second');grid on;
% figure
% for i=1:32
%     scatter(output.gpst,output.res_GPS(i,:),'.')
%     hold on
% end
% hold off
% grid on
% title('GPS residual')
% xlabel('Receiver time using GPS second');
% ylabel('Residual, unit: meter');
% figure
% scatter(output.gpst,output.sv_num_GPS,'.')
% title('Number of GPS satellites')
% xlabel('Receiver time using GPS second')
% ylabel('Distance, unit: meter');grid on
% 
% figure
% scatter(output.gpst,output.sv_num_GLO,'.')
% title('Number of GLO satellites')
% xlabel('Receiver time using GPS second')
% ylabel('Distance, unit: meter');grid on
% 
% figure
% scatter(output.gpst,output.sv_num_GAL,'.')
% title('Number of GAL satellites')
% xlabel('Receiver time using GPS second')
% ylabel('Distance, unit: meter');grid on
% 
% figure
% scatter(output.gpst,output.sv_num_BDS,'.')
% title('Number of BDS satellites')
% xlabel('Receiver time using GPS second')
% ylabel('Distance, unit: meter');grid on