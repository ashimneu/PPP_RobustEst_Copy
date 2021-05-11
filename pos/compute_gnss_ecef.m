function  log = compute_gnss_ecef(p,eph,obs)
% This function is to implement GNSS positioning with
% standard mode (without Iono, Trop, Es correction)
% or PPP mode.
% Output: log is a data struct that store the results.
%----------------------------------%
N = length(obs.tr_sow); % The number of positioning points
% Initialize output
log.gpst = obs.tr_sow-obs.tr_sow(1);
log.err = NaN(1,N); % The position (Norm) error between estimated pos and true pos
log.hor_err = NaN(1,N); % The horizontal position (Norm) error between estimated pos and true pos
log.ned_err_norm = NaN(1,N); % NED frame norm error
log.ned_err = NaN(3,N); % NED frame error
log.pos_ecef = NaN(3,N); % Estimated position in ECEF
log.v_ecef = NaN(3,N); % Estimated velocity in ECEF
log.rover_clk = NaN(1,N); % Receiver clock bias (meter)
log.clk_drift = NaN(1,N); % Receiver clock bias (meter)
log.sv_num_GPS = NaN(1,N); % The amount of GPS satellite be used
log.sv_num_GLO = NaN(1,N); % The amount of GLO satellit  used
log.sv_num_GAL = NaN(1,N); % The amount of GAL satellite be used
log.sv_num_BDS = NaN(1,N); % The amount of BDS satellite be used



if ~isempty(obs.GPS)
    log.num_obs_gps = size(obs.GPS(1).data.P,1); % The maximum of PRN recorded in obs data
else
    log.num_obs_gps = 0;
end
if ~isempty(obs.GLO)
    log.num_obs_glo = size(obs.GLO(1).data.P,1); % The maximum of PRN recorded in obs data
else
    log.num_obs_glo = 0;
end
if ~isempty(obs.GAL)
    log.num_obs_gal = size(obs.GAL(1).data.P,1); % The maximum of PRN recorded in obs data
else
    log.num_obs_gal = 0;
end
if ~isempty(obs.BDS)
    log.num_obs_bds = size(obs.BDS(1).data.P,1); % The maximum of PRN recorded in obs data
else
    log.num_obs_bds = 0;
end

%--------------------------------------------------------------------------
max_num_sv = sum([log.num_obs_gps, log.num_obs_glo,...
                log.num_obs_gal, log.num_obs_bds]);

% Initialize output
err_LS = NaN(1,N);      % position error norm
hor_err_LS = NaN(1,N);  % horizontal position error norm
ned_err_LS = NaN(1,N);  % NED frame error
GDOPLS = NaN(1,N);      % GDOP function (meter)
resLStemp = NaN(2*max_num_sv,N);
byLStemp  = NaN(2*max_num_sv,N);
H_posLStemp = NaN(2*max_num_sv,6,N);
nsvLS  = NaN(1,N); 
dnsvLS = NaN(1,N); 
npriorLS  = NaN(1,N); 
dnpriorLS = NaN(1,N);
res_stdLStemp = NaN(2*max_num_sv,N);

err_LTS = NaN(p.LTSn,N);        hor_err_LTS = NaN(p.LTSn,N);  
ned_err_LTS = NaN(p.LTSn,N);    GDOPLTS = NaN(p.LTSn,N);
resLTS = cell(p.LTSn,1);        resLTStemp = NaN(2*max_num_sv*p.LTSn,N);
nsvLTS = NaN(p.LTSn,N);         dnsvLTS = NaN(p.LTSn,N); 
npriorLTS = NaN(p.LTSn,N);      dnpriorLTS = NaN(p.LTSn,N);
byLTStemp = NaN(2*max_num_sv*p.LTSn,N);
H_posLTStemp = NaN(2*max_num_sv*p.LTSn,6,N);
res_stdLTStemp = NaN(2*max_num_sv*p.LTSn,N);

err_RAPS = NaN(p.RAPSn,N);      hor_err_RAPS = NaN(p.RAPSn,N); 
ned_err_RAPS = NaN(p.RAPSn,N);  GDOPRAPS = NaN(p.RAPSn,N);
resRAPS = cell(p.RAPSn,1);      resRAPStemp = NaN(2*max_num_sv*p.RAPSn,N);
nsvRAPS = NaN(p.RAPSn,N);       dnsvRAPS = NaN(p.RAPSn,N); 
npriorRAPS = NaN(p.RAPSn,N);    dnpriorRAPS= NaN(p.RAPSn,N);
byRAPStemp = NaN(2*max_num_sv*p.RAPSn,N);
H_posRAPStemp = NaN(2*max_num_sv*p.RAPSn,6,N);
res_stdRAPStemp = NaN(2*max_num_sv*p.RAPSn,N);

err_TD = NaN(p.TDn,N);          hor_err_TD = NaN(p.TDn,N); 
ned_err_TD = NaN(p.TDn,N);      GDOPTD = NaN(p.TDn,N);
resTD = cell(p.TDn,1);          resTDtemp = NaN(2*max_num_sv*p.TDn,N);
nsvTD = NaN(p.TDn,N);           dnsvTD = NaN(p.TDn,N); 
npriorTD = NaN(p.TDn,N);        dnpriorTD= NaN(p.TDn,N);
byTDtemp = NaN(2*max_num_sv*p.TDn,N);
H_posTDtemp = NaN(2*max_num_sv*p.TDn,6,N);
res_stdTDtemp = NaN(2*max_num_sv*p.TDn,N);

err_MShb = NaN(p.MShbn,N);      hor_err_MShb = NaN(p.MShbn,N); 
ned_err_MShb = NaN(p.MShbn,N);  GDOPMShb = NaN(p.MShbn,N);
resMShb = cell(p.MShbn,1);      resMShbtemp = NaN(2*max_num_sv*p.MShbn,N);
nsvMShb = NaN(p.MShbn,N);       dnsvMShb = NaN(p.MShbn,N); 
npriorMShb = NaN(p.MShbn,N);    dnpriorMShb= NaN(p.MShbn,N);
byMShbtemp = NaN(2*max_num_sv*p.MShbn,N);
H_posMShbtemp = NaN(2*max_num_sv*p.MShbn,6,N);
res_stdMShbtemp = NaN(2*max_num_sv*p.MShbn,N);

err_MStk = NaN(p.MStkn,N);      hor_err_MStk = NaN(p.MStkn,N); 
ned_err_MStk = NaN(p.MStkn,N);  GDOPMStk = NaN(p.MStkn,N);
resMStk = cell(p.MStkn,1);      resMStktemp = NaN(2*max_num_sv*p.MStkn,N);
nsvMStk = NaN(p.MStkn,N);       dnsvMStk = NaN(p.MStkn,N); 
npriorMStk = NaN(p.MStkn,N);    dnpriorMStk= NaN(p.MStkn,N);
byMStktemp = NaN(2*max_num_sv*p.MStkn,N);
H_posMStktemp = NaN(2*max_num_sv*p.MStkn,6,N);
res_stdMStktemp = NaN(2*max_num_sv*p.MStkn,N);

err_LSS = NaN(p.LSSn,N);        hor_err_LSS = NaN(p.LSSn,N); 
ned_err_LSS = NaN(p.LSSn,N);    GDOPLSS = NaN(p.LSSn,N);
resLSS = cell(p.LSSn,1);        resLSStemp = NaN(2*max_num_sv*p.LSSn,N);
nsvLSS = NaN(p.LSSn,N);         dnsvLSS = NaN(p.LSSn,N); 
npriorLSS = NaN(p.LSSn,N);      dnpriorLSS= NaN(p.LSSn,N);
byLSStemp = NaN(2*max_num_sv*p.LSSn,N);
H_posLSStemp = NaN(2*max_num_sv*p.LSSn,6,N);
res_stdLSStemp = NaN(2*max_num_sv*p.LSSn,N);

% For EKF -----------------------------
% state vector & covar matrix estimates at propagation steps
log.prior_state.LS       = nan(p.ns,p.numPropSteps*N); 
log.prior_statecovar.LS  = nan(p.ns,p.ns,p.numPropSteps*N);

% state vector & covar matrix estimates at measurement update steps
log.post_state.LS       = nan(p.ns,N+1); 
log.post_statecovar.LS  = nan(p.ns,p.ns,N+1);
log.delta_x.LS          = nan(p.ns,N+1);
log.Jydiag.LS           = nan(p.ns,N+1);

% initial state vector & covar matrix, i.e. at time = 0
x_post.LS = p.x_prior; % x_post.LS(1:3) = p.Grdpos.pos(1:3,1); % p.x_prior;
P_post.LS = p.P_prior;
log.post_state.LS(:,1)        = p.x_prior;
log.post_statecovar.LS(:,:,1) = p.P_prior;

if p.eb_LTS
    for idx = 1:p.LTSn
        x_post.LTS{idx} = p.x_prior;
        P_post.LTS{idx} = p.P_prior; 
        log.prior_state.LTS{idx}      = nan(p.ns,p.numPropSteps*N); 
        log.prior_statecovar.LTS{idx} = nan(p.ns,p.ns,p.numPropSteps*N);
        log.post_state.LTS{idx}       = nan(p.ns,N+1); 
        log.post_statecovar.LTS{idx}  = nan(p.ns,p.ns,N+1);
        log.post_state.LTS{idx}(:,1)        = p.x_prior;
        log.post_statecovar.LTS{idx}(:,:,1) = p.P_prior;
        log.delta_x.LTS{idx} = nan(p.ns,N+1);
        log.Jydiag.LTS{idx} = nan(p.ns,N+1);
    end
end

if p.eb_RAPS
    for idx = 1:p.RAPSn
        x_post.RAPS{idx} = p.x_prior;
        P_post.RAPS{idx} = p.P_prior; 
        log.prior_state.RAPS{idx}      = zeros(p.ns,p.numPropSteps*N); 
        log.prior_statecovar.RAPS{idx} = zeros(p.ns,p.ns,p.numPropSteps*N);
        log.post_state.RAPS{idx}       = zeros(p.ns,N+1); 
        log.post_statecovar.RAPS{idx}  = zeros(p.ns,p.ns,N+1);
        log.post_state.RAPS{idx}(:,1)        = p.x_prior;
        log.post_statecovar.RAPS{idx}(:,:,1) = p.P_prior;
        log.delta_x.RAPS{idx} = nan(p.ns,N+1);
        log.Jydiag.RAPS{idx} = nan(p.ns,N+1);
    end
end

if p.eb_TD
    for idx = 1:p.TDn
        x_post.TD{idx} = p.x_prior;
        P_post.TD{idx} = p.P_prior; 
        log.prior_state.TD{idx}      = zeros(p.ns,p.numPropSteps*N); 
        log.prior_statecovar.TD{idx} = zeros(p.ns,p.ns,p.numPropSteps*N);
        log.post_state.TD{idx}       = zeros(p.ns,N+1); 
        log.post_statecovar.TD{idx}  = zeros(p.ns,p.ns,N+1);
        log.post_state.TD{idx}(:,1)        = p.x_prior;
        log.post_statecovar.TD{idx}(:,:,1) = p.P_prior;
        log.delta_x.TD{idx} = nan(p.ns,N+1);
        log.Jydiag.TD{idx} = nan(p.ns,N+1);
    end
end

if p.eb_MShb
    for idx = 1:p.MShbn
        x_post.MShb{idx} = p.x_prior;
        P_post.MShb{idx} = p.P_prior; 
        log.prior_state.MShb{idx}      = zeros(p.ns,p.numPropSteps*N); 
        log.prior_statecovar.MShb{idx} = zeros(p.ns,p.ns,p.numPropSteps*N);
        log.post_state.MShb{idx}       = zeros(p.ns,N+1); 
        log.post_statecovar.MShb{idx}  = zeros(p.ns,p.ns,N+1);
        log.post_state.MShb{idx}(:,1)        = p.x_prior;
        log.post_statecovar.MShb{idx}(:,:,1) = p.P_prior;
        log.delta_x.MShb{idx} = nan(p.ns,N+1);
        log.Jydiag.MShb{idx} = nan(p.ns,N+1);
    end
end

if p.eb_MStk
    for idx = 1:p.MStkn
        x_post.MStk{idx} = p.x_prior;
        P_post.MStk{idx} = p.P_prior; 
        log.prior_state.MStk{idx}      = zeros(p.ns,p.numPropSteps*N); 
        log.prior_statecovar.MStk{idx} = zeros(p.ns,p.ns,p.numPropSteps*N);
        log.post_state.MStk{idx}       = zeros(p.ns,N+1); 
        log.post_statecovar.MStk{idx}  = zeros(p.ns,p.ns,N+1);
        log.post_state.MStk{idx}(:,1)        = p.x_prior;
        log.post_statecovar.MStk{idx}(:,:,1) = p.P_prior;
        log.delta_x.MStk{idx} = nan(p.ns,N+1);
        log.Jydiag.MStk{idx} = nan(p.ns,N+1);
    end
end

if p.eb_LSS
    for idx = 1:p.LSSn
        x_post.LSS{idx} = p.x_prior;
        P_post.LSS{idx} = p.P_prior; 
        log.prior_state.LSS{idx}      = zeros(p.ns,p.numPropSteps*N); 
        log.prior_statecovar.LSS{idx} = zeros(p.ns,p.ns,p.numPropSteps*N);
        log.post_state.LSS{idx}       = zeros(p.ns,N+1); 
        log.post_statecovar.LSS{idx}  = zeros(p.ns,p.ns,N+1);
        log.post_state.LSS{idx}(:,1)        = p.x_prior;
        log.post_statecovar.LSS{idx}(:,:,1) = p.P_prior;
        log.delta_x.LSS{idx} = nan(p.ns,N+1);
        log.Jydiag.LSS{idx} = nan(p.ns,N+1);
    end
end

log.xhat = nan(13,N);

log.clkbNLS = nan(1,N+1);
log.msr_count = 0;
log.N = N;
got_new_xhat_posterior = true;

log.propTime = nan(1,N*p.numPropSteps);
log.postTime = nan(1,N+1); log.postTime(1) = 0;
%--------------------------------------------------------------------------
log.msr_res_GPS = NaN(log.num_obs_gps,N); % The residual at the end
log.msr_res_GLO = NaN(log.num_obs_glo,N); % The residual at the end
log.msr_res_GAL = NaN(log.num_obs_gal,N); % The residual at the end
log.msr_res_BDS = NaN(log.num_obs_bds,N); % The residual at the end
log.msr_res = [log.msr_res_GPS;log.msr_res_GAL;log.msr_res_GLO;log.msr_res_BDS];

log.res_GPS = NaN(log.num_obs_gps,N); % The residual at the end
log.res_GLO = NaN(log.num_obs_glo,N); % The residual at the end
log.res_GAL = NaN(log.num_obs_gal,N); % The residual at the end
log.res_BDS = NaN(log.num_obs_bds,N); % The residual at the end
log.elev_GPS = NaN(log.num_obs_gps,N); % The elevation of satellites
log.elev_GLO = NaN(log.num_obs_glo,N); 
log.elev_GAL = NaN(log.num_obs_gal,N); 
log.elev_BDS = NaN(log.num_obs_bds,N); 
log.res = [log.res_GPS;log.res_GAL;log.res_GLO;log.res_BDS];
log.elev = [log.elev_GPS;log.elev_GAL;log.elev_GLO;log.elev_BDS];
log.t_kf = []; % Time propagation stamp
log.state_PVA = []; % State data
log.P_cov = []; % Covariance data
% Mark the sat prn that be computed
gpslog.svprn_mark = zeros(log.num_obs_gps,1);glolog.svprn_mark = zeros(log.num_obs_glo,1);
gallog.svprn_mark = zeros(log.num_obs_gal,1);bdslog.svprn_mark = zeros(log.num_obs_bds,1);
% Record the prn of each system for whose satllite been used
gpslog.prn_record = zeros(log.num_obs_gps,1);glolog.prn_record = zeros(log.num_obs_glo,1);
gallog.prn_record = zeros(log.num_obs_gal,1);bdslog.prn_record = zeros(log.num_obs_bds,1);
% Satellite position in ECEF frame
gpslog.s_pos_ecef = [];glolog.s_pos_ecef = [];gallog.s_pos_ecef = [];bdslog.s_pos_ecef = [];
if p.post_mode == 1
% precise satellite position in ECEF frame
gpslog.s_pos_prc = [];glolog.s_pos_prc = [];gallog.s_pos_prc = [];bdslog.s_pos_prc = [];
end
% Satellite velocity in ECEF frame
gpslog.s_v_ecef = [];glolog.s_v_ecef = [];gallog.s_v_ecef = [];bdslog.s_v_ecef = [];
% Signal propagation time
gpslog.tp = [];glolog.tp = [];gallog.tp = [];bdslog.tp = [];
% corrected pseudorange
gpslog.corr_range = [];glolog.corr_range = [];gallog.corr_range = [];bdslog.corr_range = [];
% Doppler
gpslog.dp_range = [];glolog.dp_range = [];gallog.dp_range = [];bdslog.dp_range = [];
% The number of satellite be computed
gpslog.num_sv = 0;glolog.num_sv = 0;gallog.num_sv = 0;bdslog.num_sv = 0;

log.NN = N;
for i = 1:p.inval:log.NN
%--------------------------------------------------------------------------
    sidx = (i-1)*p.numPropSteps + 1;  % start index
    eidx = i*p.numPropSteps;          % end index
    
    log.propTime(sidx:eidx) = (i-1) + p.T : p.T : i;
    % batch - bundle of state vector between successive measurement updates
    % [Statebatch(i).xLS,Covarbatch(i).PLS]= KFpropagate(p,x_post,P_post); % all batches of state
    if got_new_xhat_posterior
        [log.prior_state.LS(:,sidx:eidx),log.prior_statecovar.LS(:,:,sidx:eidx)]= KFpropagate(p,x_post.LS,P_post.LS);
    else
        [log.prior_state.LS(:,sidx:eidx),log.prior_statecovar.LS(:,:,sidx:eidx)]= KFpropagate(p,x_prior.LS,P_prior.LS);
    end
    
    % state & Covariance after some propagations
    x_prior.LS = log.prior_state.LS(:,eidx); 
    P_prior.LS = log.prior_statecovar.LS(:,:,eidx);
    
    if p.eb_LTS
    for idx = 1:p.LTSn
        if got_new_xhat_posterior
            [log.prior_state.LTS{idx}(:,sidx:eidx),log.prior_statecovar.LTS{idx}(:,:,sidx:eidx)]= KFpropagate(p,x_post.LTS{idx},P_post.LTS{idx});
        else
            [log.prior_state.LTS{idx}(:,sidx:eidx),log.prior_statecovar.LTS{idx}(:,:,sidx:eidx)]= KFpropagate(p,x_prior.LTS{idx},P_prior.LTS{idx});
        end
        x_prior.LTS{idx} = log.prior_state.LTS{idx}(:,eidx); 
        P_prior.LTS{idx} = log.prior_statecovar.LTS{idx}(:,:,eidx);
    end
    end
    
    if p.eb_RAPS
    for idx = 1:p.RAPSn
        if got_new_xhat_posterior
            [log.prior_state.RAPS{idx}(:,sidx:eidx),log.prior_statecovar.RAPS{idx}(:,:,sidx:eidx)]= KFpropagate(p,x_post.RAPS{idx},P_post.RAPS{idx});
        else
            [log.prior_state.RAPS{idx}(:,sidx:eidx),log.prior_statecovar.RAPS{idx}(:,:,sidx:eidx)]= KFpropagate(p,x_prior.RAPS{idx},P_prior.RAPS{idx});
        end
        x_prior.RAPS{idx} = log.prior_state.RAPS{idx}(:,eidx); 
        P_prior.RAPS{idx} = log.prior_statecovar.RAPS{idx}(:,:,eidx);
    end
    end
    
    if p.eb_TD
    for idx = 1:p.TDn
        if got_new_xhat_posterior
            [log.prior_state.TD{idx}(:,sidx:eidx),log.prior_statecovar.TD{idx}(:,:,sidx:eidx)]= KFpropagate(p,x_post.TD{idx},P_post.TD{idx});
        else
            [log.prior_state.TD{idx}(:,sidx:eidx),log.prior_statecovar.TD{idx}(:,:,sidx:eidx)]= KFpropagate(p,x_prior.TD{idx},P_prior.TD{idx});
        end
        x_prior.TD{idx} = log.prior_state.TD{idx}(:,eidx); 
        P_prior.TD{idx} = log.prior_statecovar.TD{idx}(:,:,eidx);
    end
    end
    
    if p.eb_MShb
    for idx = 1:p.MShbn
        if got_new_xhat_posterior
            [log.prior_state.MShb{idx}(:,sidx:eidx),log.prior_statecovar.MShb{idx}(:,:,sidx:eidx)]= KFpropagate(p,x_post.MShb{idx},P_post.MShb{idx});
        else
            [log.prior_state.MShb{idx}(:,sidx:eidx),log.prior_statecovar.MShb{idx}(:,:,sidx:eidx)]= KFpropagate(p,x_prior.MShb{idx},P_prior.MShb{idx});
        end
        x_prior.MShb{idx} = log.prior_state.MShb{idx}(:,eidx); 
        P_prior.MShb{idx} = log.prior_statecovar.MShb{idx}(:,:,eidx);
    end
    end
    
    if p.eb_MStk
    for idx = 1:p.MStkn
        if got_new_xhat_posterior
            [log.prior_state.MStk{idx}(:,sidx:eidx),log.prior_statecovar.MStk{idx}(:,:,sidx:eidx)]= KFpropagate(p,x_post.MStk{idx},P_post.MStk{idx});
        else
            [log.prior_state.MStk{idx}(:,sidx:eidx),log.prior_statecovar.MStk{idx}(:,:,sidx:eidx)]= KFpropagate(p,x_prior.MStk{idx},P_prior.MStk{idx});
        end
        x_prior.MStk{idx} = log.prior_state.MStk{idx}(:,eidx); 
        P_prior.MStk{idx} = log.prior_statecovar.MStk{idx}(:,:,eidx);
    end
    end
    
    if p.eb_LSS
    for idx = 1:p.LSSn
        if got_new_xhat_posterior
            [log.prior_state.LSS{idx}(:,sidx:eidx),log.prior_statecovar.LSS{idx}(:,:,sidx:eidx)]= KFpropagate(p,x_post.LSS{idx},P_post.LSS{idx});
        else
            [log.prior_state.LSS{idx}(:,sidx:eidx),log.prior_statecovar.LSS{idx}(:,:,sidx:eidx)]= KFpropagate(p,x_prior.LSS{idx},P_prior.LSS{idx});
        end
        x_prior.LSS{idx} = log.prior_state.LSS{idx}(:,eidx); 
        P_prior.LSS{idx} = log.prior_statecovar.LSS{idx}(:,:,eidx);
    end
    end
    
%--------------------------------------------------------------------------
    if ~isnan(p.Grdpos.t(1))
        index = round(p.Grdpos.t) == round(obs.tr_sow(i));
        grdpos = p.Grdpos.pos(:,index);
        if isempty(grdpos)
            continue;
        end
    else
        grdpos = p.Grdpos.pos;
    end
    if mod(i,4000)==0
        i
    end
    
    % Time propagation
    if i>1 && p.PVA_enable == 1
        t = p.t(i-1):seconds(0.1):p.t(i);
        for k = 1:length(t)-1
            dt = seconds(t(k+1)-t(k));
            log.t_kf = [log.t_kf,t(k+1)];
            p = time_prop(p,dt);
            log.state_PVA = [log.state_PVA,p.state_PVA];
            log.P_cov = [log.P_cov,diag(p.P_cov)];
        end
    end
    
    p.i = i; % To debug
    if ~isempty(obs.GPS)&& p.freq==1 && p.enableGPS == 1
        % GPS satellite position computation, Single frenquency receiver mode    
        gpslog = satpost_corrpsedR_singlefreq(p,eph,obs,i,log.num_obs_gps,'GPS');
    end
    if ~isempty(obs.GLO)&& p.freq==1 && p.enableGLO == 1
        % GLO satellite position computation, Single frenquency receiver mode
        glolog = satpost_corrpsedR_singlefreq(p,eph,obs,i,log.num_obs_glo,'GLO');
    end
    if ~isempty(obs.GAL)&& p.freq==1 && p.enableGAL == 1
        % GAL satellite position computation, Single frenquency receiver mode
        gallog = satpost_corrpsedR_singlefreq(p,eph,obs,i,log.num_obs_gal,'GAL');
    end
    if ~isempty(obs.BDS)&& p.freq==1 && p.enableBDS == 1
        % BDS satellite position computation, Single frenquency receiver mode
        bdslog = satpost_corrpsedR_singlefreq(p,eph,obs,i,log.num_obs_bds,'BDS');
    end
    cpt.prn_record = [gpslog.prn_record;glolog.prn_record;gallog.prn_record;bdslog.prn_record];
    cpt.svprn_mark = [gpslog.svprn_mark;glolog.svprn_mark;gallog.svprn_mark;bdslog.svprn_mark];
    cpt.corr_range = [gpslog.corr_range;glolog.corr_range;gallog.corr_range;bdslog.corr_range];
    cpt.dp_range = [gpslog.dp_range;glolog.dp_range;gallog.dp_range;bdslog.dp_range];
    ind = find(cpt.corr_range==0);
    cpt.corr_range(ind) = [];cpt.dp_range(ind) = [];
    cpt.s_pos_ecef = [gpslog.s_pos_ecef,glolog.s_pos_ecef,gallog.s_pos_ecef,bdslog.s_pos_ecef];
    cpt.s_pos_ecef(:,ind) = [];
    if p.post_mode == 1
        cpt.s_pos_prc = [gpslog.s_pos_prc,glolog.s_pos_prc,gallog.s_pos_prc,bdslog.s_pos_prc];
        cpt.s_pos_prc(:,ind) = [];
    end
    cpt.s_v_ecef = [gpslog.s_v_ecef,glolog.s_v_ecef,gallog.s_v_ecef,bdslog.s_v_ecef];
    cpt.s_v_ecef(:,ind) = [];
    cpt.tp = [gpslog.tp;glolog.tp;gallog.tp;bdslog.tp];
    cpt.tp(ind) = [];
    cpt.num_sv = [gpslog.num_sv,glolog.num_sv,gallog.num_sv,bdslog.num_sv];
    % elevation & azimuth
    cpt.elev = NaN(length(cpt.corr_range),1); cpt.az = NaN(length(cpt.corr_range),1);
    % trop delay and iono delay
    cpt.trop_delay = NaN(length(cpt.corr_range),1); cpt.iono_delay = NaN(length(cpt.corr_range),1);
    
    if sum(cpt.num_sv)>=p.min_sv
        if p.state_PVA(1)==0
            for kk = 1:3
                [p,~] = userpos(p,cpt);
                if p.post_mode == 1
                    p.mk = 1;
                end
            end
            log.t_kf = [log.t_kf,p.t(i)];
            log.state_PVA = [log.state_PVA,p.state_PVA];
            log.P_cov = [log.P_cov,diag(p.P_cov)];
            
        end
        % Rotate the sat pos to common reference frame
        cpt = earth_rotation_corr(p,cpt,p.state_PVA(10)/p.c);
        % Check elevation
        cpt = elevaz_check(p,cpt,p.state_PVA(1:3));
        if sum(cpt.num_sv)>=p.min_sv
            switch p.post_mode
                case 0 % Standard GNSS
%                     if p.elev_mark ==0
%                         log = save_result(p,cpt,log,i,re_pos,clock_bias,res);
%                     else
                        % Compute the final position
%                       [re_pos,clock_bias,res] = userpos_Rcorr(p,cpt);
                      tdoy = doy(obs.tr_prime(1:3,i)); % Day of year
                      [rt.week, rt.dow, rt.sow] = date2gnsst(obs.tr_prime(:,i)');
                      cpt.IoFac = zeros(length(cpt.corr_range),1);
                      cpt = trop_iono_compute(p,eph,cpt,obs,p.state0(1:3),tdoy,[],rt);
                      cpt.corr_range = cpt.corr_range - cpt.trop_delay - cpt.iono_delay;
                      [re_pos,clock_bias,res] = userpos(p,cpt);
                      [log,p.state0] = save_result(p,cpt,log,i,re_pos,clock_bias,res,grdpos);
%                     end
                case 1 % PPP
                    tdoy = doy(obs.tr_prime(1:3,i)); % Day of year
                    [rt.week, rt.dow, rt.sow] = date2gnsst(obs.tr_prime(:,i)');
                    %%-------------%%
                    cpt.IoFac = zeros(length(cpt.corr_range),1);
                    cpt = trop_iono_compute(p,eph,cpt,obs,p.state_PVA(1:3),tdoy,p.USTEC,rt);
                    % Using the correction to the measurements
                    if ~isempty(find(cpt.iono_delay~=0, 1))
                    cpt.corr_range = cpt.corr_range - cpt.trop_delay - cpt.iono_delay;
                    %%-------------%%
                    % Compute the final position
%                         [p,res] = userpos(p,cpt);
                    p.state_PVA_1 = p.state_PVA;
                    p.P_cov_1 = p.P_cov;
                    
                    if i==1
                        [p,res] = userpos(p,cpt);
                        p.state_PVA_1 = p.state_PVA;
                        p.P_cov_1 = p.P_cov;
                    else
                        [p,res] = EKF_update(p,cpt);
%                             [p,res] = userpos(p,cpt);
                    end
                    log.state_PVA(:,end) = p.state_PVA;
                    log.P_cov(:,end) = diag(p.P_cov);
                    log = save_result(p,cpt,log,i,res,grdpos);
                    
                    
                    end
%--------------------------------------------------------------------------
                    ind_mark = cpt.svprn_mark ~= 0;
                    
%                     if p.i == 1
%                         x_prior.LS = p.state_PVA_1; 
%                         P_prior.LS = p.P_cov_1; 
%                         x_prior.LTS{1} = p.state_PVA_1; 
%                         P_prior.LTS{1} = p.P_cov_1;
%                     end
                    
                    [re_pos,msr_res,cpt.GDOP,cpt.nsv,cpt.dnsv,cpt.nprior,cpt.dnprior,postxhat,...
                    postcovar,delta_x,Jydiag,by,H_pos,cpt.outliervec,cpt.outlierbin,cpt.res_std] = ...
                    userposlinear_PVA(p,cpt,grdpos,x_prior,P_prior);
                    log.msr_count = log.msr_count + 1; % number of measurement updates
                    got_new_xhat_posterior = true;
                    log.postTime(i+1) = i; % GPS time when measurement update process is executed.
                    %------------Save Generated Outliers to Database------%
                    log.outliervec{i} = cpt.outliervec;
                    log.outlierbin{i} = cpt.outlierbin;
                    %-----------------------PVA---------------------------%
                    x_post.LS = postxhat.LS;   
                    P_post.LS = postcovar.LS;
                    log.post_state.LS(:,i+1)        = postxhat.LS;
                    log.post_statecovar.LS(:,:,i+1) = postcovar.LS;
                    log.delta_x.LS(:,i+1)           = delta_x.LS;
                    log.Jydiag.LS(:,i+1)            = Jydiag.LS;
                    
                    if p.eb_LTS == 1
                        for idx = 1:p.LTSn
                            x_post.LTS{idx} = postxhat.LTS{idx};
                            P_post.LTS{idx} = postcovar.LTS{idx};
                            log.post_state.LTS{idx}(:,i+1) = postxhat.LTS{idx};
                            log.post_statecovar.LTS{idx}(:,:,i+1) = postcovar.LTS{idx};
                            log.delta_x.LTS{idx}(:,i+1) = delta_x.LTS{idx};
                            log.Jydiag.LTS{idx}(:,i+1)  = Jydiag.LTS{idx};
                        end
                    end
                    
                    if p.eb_RAPS == 1
                        for idx = 1:p.RAPSn
                            x_post.RAPS{idx} = postxhat.RAPS{idx};
                            P_post.RAPS{idx} = postcovar.RAPS{idx};
                            log.post_state.RAPS{idx}(:,i+1) = postxhat.RAPS{idx};
                            log.post_statecovar.RAPS{idx}(:,:,i+1) = postcovar.RAPS{idx};
                            log.delta_x.RAPS{idx}(:,i+1) = delta_x.RAPS{idx};
                            log.Jydiag.RAPS{idx}(:,i+1)  = Jydiag.RAPS{idx};
                        end
                    end
                    
                    if p.eb_TD == 1
                        for idx = 1:p.TDn
                            x_post.TD{idx} = postxhat.TD{idx};
                            P_post.TD{idx} = postcovar.TD{idx};
                            log.post_state.TD{idx}(:,i+1) = postxhat.TD{idx};
                            log.post_statecovar.TD{idx}(:,:,i+1) = postcovar.TD{idx};
                            log.delta_x.TD{idx}(:,i+1) = delta_x.TD{idx};
                            log.Jydiag.TD{idx}(:,i+1)  = Jydiag.TD{idx};
                        end
                    end
                    
                    if p.eb_MShb == 1
                        for idx = 1:p.MShbn
                            x_post.MShb{idx} = postxhat.MShb{idx};
                            P_post.MShb{idx} = postcovar.MShb{idx};
                            log.post_state.MShb{idx}(:,i+1) = postxhat.MShb{idx};
                            log.post_statecovar.MShb{idx}(:,:,i+1) = postcovar.MShb{idx};
                            log.delta_x.MShb{idx}(:,i+1) = delta_x.MShb{idx};
                            log.Jydiag.MShb{idx}(:,i+1)  = Jydiag.MShb{idx};
                        end
                    end
                    
                    if p.eb_MStk == 1
                        for idx = 1:p.MStkn
                            x_post.MStk{idx} = postxhat.MStk{idx};
                            P_post.MStk{idx} = postcovar.MStk{idx};
                            log.post_state.MStk{idx}(:,i+1) = postxhat.MStk{idx};
                            log.post_statecovar.MStk{idx}(:,:,i+1) = postcovar.MStk{idx};
                            log.delta_x.MStk{idx}(:,i+1) = delta_x.MStk{idx};
                            log.Jydiag.MStk{idx}(:,i+1)  = Jydiag.MStk{idx};
                        end
                    end
                    
                    if p.eb_LSS == 1
                        for idx = 1:p.LSSn
                            x_post.LSS{idx} = postxhat.LSS{idx};
                            P_post.LSS{idx} = postcovar.LSS{idx};
                            log.post_state.LSS{idx}(:,i+1) = postxhat.LSS{idx};
                            log.post_statecovar.LSS{idx}(:,:,i+1) = postcovar.LSS{idx};
                            log.delta_x.LSS{idx}(:,i+1) = delta_x.LSS{idx};
                            log.Jydiag.LSS{idx}(:,i+1)  = Jydiag.LSS{idx};
                        end
                    end
%--------------------------------------------------------------------------                    
                    [ned_err_LS(:,i),hor_err_LS(:,i),err_LS(:,i),GDOPLS(:,i),...
                    nsvLS(i),dnsvLS(i),npriorLS(i),dnpriorLS(i),resLStemp(:,i),...
                    byLStemp(:,i),H_posLStemp(:,:,i),res_stdLStemp(:,i)] = ...
                    save_errNorm_res_GDOP(1,grdpos,re_pos.LS,cpt.GDOP.LS,...
                    cpt.nsv.LS,cpt.dnsv.LS,cpt.nprior.LS,cpt.dnprior.LS,...
                    max_num_sv,ind_mark,msr_res.LS,by.LS,H_pos.LS,cpt.res_std.LS);

                                      
                    if p.eb_LTS == 1
                        [ned_err_LTS(:,i),hor_err_LTS(:,i),err_LTS(:,i),GDOPLTS(:,i),...
                        nsvLTS(:,i),dnsvLTS(:,i),npriorLTS(:,i),dnpriorLTS(:,i),...
                        resLTStemp(:,i),byLTStemp(:,i),H_posLTStemp(:,:,i),res_stdLTStemp(:,i)] = save_errNorm_res_GDOP(p.LTSn,...
                        grdpos,re_pos.LTS,cpt.GDOP.LTS,cpt.nsv.LTS,cpt.dnsv.LTS,...
                        cpt.nprior.LTS,cpt.dnprior.LTS,max_num_sv,ind_mark,msr_res.LTS,by.LTS,H_pos.LTS,cpt.res_std.LTS);                                                
                    end
                    
                    if p.eb_RAPS == 1 
                        [ned_err_RAPS(:,i),hor_err_RAPS(:,i),err_RAPS(:,i),GDOPRAPS(:,i),...
                        nsvRAPS(:,i),dnsvRAPS(:,i),npriorRAPS(:,i),dnpriorRAPS(:,i),...
                        resRAPStemp(:,i),byRAPStemp(:,i),H_posRAPStemp(:,:,i),res_stdRAPStemp(:,i)] = save_errNorm_res_GDOP(p.RAPSn,...
                        grdpos,re_pos.RAPS,cpt.GDOP.RAPS,cpt.nsv.RAPS,cpt.dnsv.RAPS,...
                        cpt.nprior.RAPS,cpt.dnprior.RAPS,max_num_sv,ind_mark,msr_res.RAPS,by.RAPS,H_pos.RAPS,cpt.res_std.RAPS);        
                    end
                    
                    if p.eb_TD == 1
                        [ned_err_TD(:,i),hor_err_TD(:,i),err_TD(:,i),GDOPTD(:,i),...
                        nsvTD(:,i),dnsvTD(:,i),npriorTD(:,i),dnpriorTD(:,i),...
                        resTDtemp(:,i),byTDtemp(:,i),H_posTDtemp(:,:,i),res_stdTDtemp(:,i)] = save_errNorm_res_GDOP(p.TDn,...
                        grdpos,re_pos.TD,cpt.GDOP.TD,cpt.nsv.TD,cpt.dnsv.TD,...
                        cpt.nprior.TD,cpt.dnprior.TD,max_num_sv,ind_mark,msr_res.TD,by.TD,H_pos.TD,cpt.res_std.TD);                                           
                    end
                    
                    if p.eb_MShb == 1                        
                        [ned_err_MShb(:,i),hor_err_MShb(:,i),err_MShb(:,i),GDOPMShb(:,i),...
                        nsvMShb(:,i),dnsvMShb(:,i),npriorMShb(:,i),dnpriorMShb(:,i),...
                        resMShbtemp(:,i),byMShbtemp(:,i),H_posMShbtemp(:,:,i),res_stdMShbtemp(:,i)] = ...
                        save_errNorm_res_GDOP(p.MShbn,grdpos,re_pos.MShb,...
                        cpt.GDOP.MShb,cpt.nsv.MShb,cpt.dnsv.MShb,cpt.nprior.MShb,...
                        cpt.dnprior.MShb,max_num_sv,ind_mark,msr_res.MShb,by.MShb,H_pos.MShb,cpt.res_std.MShb);                                            
                    end
                    
                    if p.eb_MStk == 1
                        [ned_err_MStk(:,i),hor_err_MStk(:,i),err_MStk(:,i),GDOPMStk(:,i),...
                        nsvMStk(:,i),dnsvMStk(:,i),npriorMStk(:,i),dnpriorMStk(:,i),...
                        resMStktemp(:,i),byMStktemp(:,i),H_posMStktemp(:,:,i),res_stdMStktemp(:,i)] = save_errNorm_res_GDOP(p.MStkn,...
                        grdpos,re_pos.MStk,cpt.GDOP.MStk,cpt.nsv.MStk,cpt.dnsv.MStk,...
                        cpt.nprior.MStk,cpt.dnprior.MStk,max_num_sv,ind_mark,msr_res.MStk,by.MStk,H_pos.MStk,cpt.res_std.MStk);                        
                    end
                    
                    if p.eb_LSS == 1
                        [ned_err_LSS(:,i),hor_err_LSS(:,i),err_LSS(:,i),GDOPLSS(:,i),...
                        nsvLSS(:,i),dnsvLSS(:,i),npriorLSS(:,i),dnpriorLSS(:,i),...
                        resLSStemp(:,i),byLSStemp(:,i),H_posLSStemp(:,:,i),res_stdLSStemp(:,i)] =...
                        save_errNorm_res_GDOP(p.LSSn,grdpos,re_pos.LSS,...
                        cpt.GDOP.LSS,cpt.nsv.LSS,cpt.dnsv.LSS,cpt.nprior.LSS,...
                        cpt.dnprior.LSS,max_num_sv,ind_mark,msr_res.LSS,by.LSS,H_pos.LSS,cpt.res_std.LSS);
                    end

%--------------------------------------------------------------------------                        
                case 2 % DGNSS
                    if ~isempty(p.eph_b) && ~isempty(p.obs_b)
                        [cpt,n] = diff_corr_compute(p,cpt,obs.tr_sow(i));
                        if ~isempty(n) && sum(cpt.num_sv)>=p.min_sv
                            cpt = cpt_clear(cpt); % Clear the data where don't have diff correction
                            cpt.corr_range = cpt.corr_range - cpt.diff_corr;
%                             [re_pos,clock_bias,res] = userpos(p,cpt.sat_pos_Rcorr,...
%                                 cpt.corr_range,cpt.num_sv);
                        if i==1
                            [p,res] = userpos(p,cpt);
                        else
                            [p,res] = EFK_update(p,cpt);
%                             [p,res] = userpos(p,cpt);
                        end
                        log.state_PVA(:,end) = p.state_PVA;
                        log.P_cov(:,end) = diag(p.P_cov);
%                             if ~isempty(re_pos)
                            log = save_result(p,cpt,log,i,res,grdpos);
%                             end
                        end
                    else
                        warning('No differential source given')
                    end
                    
                case 3 % VRS
                    tdoy = doy(obs.tr_prime(1:3,i)); % Day of year
                    [rt.week, rt.dow, rt.sow] = date2gnsst(obs.tr_prime(:,i)');
                    rt.sow = round(rt.sow);
                    cpt = vrs_corr_compute(p,cpt,eph,tdoy,rt);
                    if length(cpt.corr_range)>=p.min_sv
                        cpt.corr_range = cpt.corr_range - cpt.diff_corr;
                        [re_pos,clock_bias,res] = userpos(p,cpt);
                        if ~isempty(re_pos)
                            [log,p.state0] = save_result(p,cpt,log,i,re_pos,clock_bias,res,grdpos);
                        end
                    end
                otherwise
                    warning('Unsupport positioning option');
            end
       end

    end
    
end

if p.genOutlier == 1 && p.saveOutlier == 1
    saveoutlier(p.outlierdbpath,log.outliervec,log.outlierbin,p.outlierparam);
end

log.err_LS  = err_LS; log.hor_err_LS = hor_err_LS;
log.ned_err_LS = ned_err_LS;
log.GDOPLS = GDOPLS;
log.nsvLS = nsvLS;
log.dnsvLS = dnsvLS;
log.npriorLS = npriorLS;
log.dnpriorLS = dnpriorLS;
log.resLS = {resLStemp};
log.res_stdLS = {res_stdLStemp};
log.byLS = {byLStemp};
log.H_posLS = {H_posLStemp};

if p.eb_RAPS == 1
    pointer1 = 1;
    for idx = 1:p.RAPSn
        resRAPS{idx} = resRAPStemp(pointer1:pointer1+2*max_num_sv-1,:);
        byRAPS{idx} = byRAPStemp(pointer1:pointer1+2*max_num_sv-1,:);
        H_posRAPS{idx} = H_posRAPStemp(pointer1:pointer1+2*max_num_sv-1,:,:);
        res_stdRAPS{idx} = res_stdRAPStemp(pointer1:pointer1+2*max_num_sv-1,:);
        pointer1 = pointer1 + 2*max_num_sv;
    end
    log.err_RAPS = err_RAPS; log.hor_err_RAPS = hor_err_RAPS;
    log.ned_err_RAPS = ned_err_RAPS;
    log.GDOPRAPS = GDOPRAPS;
    log.nsvRAPS  = nsvRAPS;
    log.dnsvRAPS = dnsvRAPS;
    log.npriorRAPS = npriorRAPS;
    log.dnpriorRAPS = dnpriorRAPS;
    log.resRAPS = resRAPS;
    log.res_stdRAPS = res_stdRAPS;
    log.byRAPS  = byRAPS;
    log.H_posRAPS = H_posRAPS;
    log.RAPSEps = p.RAPSEps; log.RAPSn = p.RAPSn;
end


if p.eb_LTS == 1
    pointer1 = 1;
    for idx = 1:p.LTSn
        resLTS{idx} = resLTStemp(pointer1:pointer1+2*max_num_sv-1,:);
        byLTS{idx}  = byLTStemp(pointer1:pointer1+2*max_num_sv-1,:);
        H_posLTS{idx} = H_posLTStemp(pointer1:pointer1+2*max_num_sv-1,:,:);
        res_stdLTS{idx} = res_stdLTStemp(pointer1:pointer1+2*max_num_sv-1,:);
        pointer1 = pointer1 + 2*max_num_sv;
    end
    log.err_LTS = err_LTS; log.hor_err_LTS = hor_err_LTS; 
    log.ned_err_LTS = ned_err_LTS;
    log.GDOPLTS = GDOPLTS;
    log.nsvLTS  = nsvLTS;
    log.dnsvLTS = dnsvLTS;
    log.npriorLTS = npriorLTS;
    log.dnpriorLTS = dnpriorLTS;
    log.resLTS = resLTS;
    log.res_stdLTS = res_stdLTS;
    log.byLTS = byLTS;
    log.H_posLTS = H_posLTS;
    log.LTSOption = p.LTSOption; log.LTSn = p.LTSn;
end

if p.eb_TD == 1
    pointer1 = 1;
    for idx = 1:p.TDn
        resTD{idx} = resTDtemp(pointer1:pointer1+2*max_num_sv-1,:);
        byTD{idx}  = byTDtemp(pointer1:pointer1+2*max_num_sv-1,:);
        H_posTD{idx} = H_posTDtemp(pointer1:pointer1+2*max_num_sv-1,:,:);
        res_stdTD{idx} = res_stdTDtemp(pointer1:pointer1+2*max_num_sv-1,:);
        pointer1 = pointer1 + 2*max_num_sv;
    end
    log.err_TD = err_TD; log.hor_err_TD = hor_err_TD; 
    log.ned_err_TD = ned_err_TD;
    log.GDOPTD = GDOPTD;
    log.nsvTD  = nsvTD;
    log.dnsvTD = dnsvTD;
    log.npriorTD = npriorTD;
    log.dnpriorTD = dnpriorTD;
    log.resTD = resTD;
    log.res_stdTD = res_stdTD;
    log.byTD = byTD;
    log.H_posTD = H_posTD;
    log.TDLambda = p.TDLambda; log.TDn = p.TDn;
end

if p.eb_MShb == 1
    pointer1 = 1;
    for idx = 1:p.MShbn
        resMShb{idx} = resMShbtemp(pointer1:pointer1+2*max_num_sv-1,:);
        byMShb{idx}  = byMShbtemp(pointer1:pointer1+2*max_num_sv-1,:);
        H_posMShb{idx} = H_posMShbtemp(pointer1:pointer1+2*max_num_sv-1,:,:);
        res_stdMShb{idx} = res_stdMShbtemp(pointer1:pointer1+2*max_num_sv-1,:);
        pointer1 = pointer1 + 2*max_num_sv;
    end
    log.err_MShb = err_MShb; log.hor_err_MShb = hor_err_MShb; 
    log.ned_err_MShb = ned_err_MShb;
    log.GDOPMShb = GDOPMShb;
    log.nsvMShb  = nsvMShb;
    log.dnsvMShb = dnsvMShb;
    log.npriorMShb = npriorMShb;
    log.dnpriorMShb = dnpriorMShb;
    log.resMShb = resMShb;
    log.res_stdMShb = res_stdMShb;
    log.byMShb = byMShb;
    log.H_posMShb = H_posMShb;
    log.MShbConst = p.MShbConst; log.MShbn = p.MShbn;
end

if p.eb_MStk == 1
    pointer1 = 1;
    for idx = 1:p.MStkn
        resMStk{idx} = resMStktemp(pointer1:pointer1+2*max_num_sv-1,:);
        byMStk{idx}  = byMStktemp(pointer1:pointer1+2*max_num_sv-1,:);
        H_posMStk{idx} = H_posMStktemp(pointer1:pointer1+2*max_num_sv-1,:,:);
        res_stdMStk{idx} = res_stdMStktemp(pointer1:pointer1+2*max_num_sv-1,:);
        pointer1 = pointer1 + 2*max_num_sv;
    end
    log.err_MStk = err_MStk; log.hor_err_MStk = hor_err_MStk; 
    log.ned_err_MStk = ned_err_MStk;
    log.GDOPMStk = GDOPMStk;
    log.nsvMStk  = nsvMStk;
    log.dnsvMStk = dnsvMStk;
    log.npriorMStk = npriorMStk;
    log.dnpriorMStk = dnpriorMStk;
    log.resMStk = resMStk;
    log.res_stdMStk = res_stdMStk;
    log.byMStk = byMStk;
    log.H_posMStk = H_posMStk;
    log.MStkConst = p.MStkConst; log.MStkn = p.MStkn;
end

if p.eb_LSS == 1
    pointer1 = 1;
    for idx = 1:p.LSSn
        resLSS{idx} = resLSStemp(pointer1:pointer1+2*max_num_sv-1,:);
        byLSS{idx}  = byLSStemp(pointer1:pointer1+2*max_num_sv-1,:);
        H_posLSS{idx} = H_posLSStemp(pointer1:pointer1+2*max_num_sv-1,:,:);
        res_stdLSS{idx} = res_stdLSStemp(pointer1:pointer1+2*max_num_sv-1,:);
        pointer1 = pointer1 + 2*max_num_sv;
    end
    log.err_LSS = err_LSS; log.hor_err_LSS = hor_err_LSS; 
    log.ned_err_LSS = ned_err_LSS;
    log.GDOPLSS = GDOPLSS;
    log.nsvLSS  = nsvLSS;
    log.dnsvLSS = dnsvLSS;
    log.npriorLSS = npriorLSS;
    log.dnpriorLSS = dnpriorLSS;
    log.resLSS = resLSS;
    log.byLSS = byLSS;
    log.H_posLSS = H_posLSS;
    log.res_stdLSS = res_stdLSS;
    log.LSSLambda = p.LSSLambda; 
    log.LSSn = p.LSSn;
end

log.p = p;
end