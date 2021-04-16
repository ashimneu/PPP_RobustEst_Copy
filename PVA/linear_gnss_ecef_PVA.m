function  log = linear_gnss_ecef_PVA(p,eph,obs)
% This function is to implement GNSS positioning with
% linear mode (without Iono, Trop, Es correction)
% Output: log is a data struct that store the results.
%----------------------------------%
N = length(obs.tr_sow);         % The number of positioning points
p.MStkn = numel(p.MStkConst); 
p.MShbn = numel(p.MShbConst); 
p.LTSn  = numel(p.LTSOption);
p.TDn   = numel(p.TDLambda);
p.RAPSn = size(p.RAPSEps,1);
p.LSSn  = numel(p.LSSLambda);
pos_ecef   = cell(1,N); % Estimated position in ECEF
rover_clk  = cell(1,N); % Receiver clock bias (meter)
sv_num_GPS = NaN(1,N);  % The amount of GPS satellite be used
sv_num_GLO = NaN(1,N);  % The amount of GLO satellite be used
sv_num_GAL = NaN(1,N);  % The amount of GAL satellite be used
sv_num_BDS = NaN(1,N);  % The amount of BDS satellite be used
%--------------------------------------------------------------------------
if ~isempty(obs.GPS)
    num_obs_gps = size(obs.GPS.P1,1); % The maximum of PRN recorded in obs data
else
    num_obs_gps = 0;
end
if ~isempty(obs.GLO)
    num_obs_glo = size(obs.GLO.P1,1); % The maximum of PRN recorded in obs data
else
    num_obs_glo = 0;
end
if ~isempty(obs.GAL)
    num_obs_gal = size(obs.GAL.P1,1); % The maximum of PRN recorded in obs data
else
    num_obs_gal = 0;
end
if ~isempty(obs.BDS)
    num_obs_bds = size(obs.BDS.P1,1); % The maximum of PRN recorded in obs data
else
    num_obs_bds = 0;
end
res_GPS = NaN(num_obs_gps,N); % The residual at the end
res_GLO = NaN(num_obs_glo,N); % The residual at the end
res_GAL = NaN(num_obs_gal,N); % The residual at the end
res_BDS = NaN(num_obs_bds,N); % The residual at the end
max_num_sv = size([res_GPS; res_GAL; res_GLO;res_BDS],1);

% Initialize output
err_LS = NaN(1,N);      % position error norm
hor_err_LS = NaN(1,N);  % horizontal position error norm
ned_err_LS = NaN(1,N);  % NED frame error
GDOPLS = NaN(1,N);      % GDOP function (meter)
resLStemp = NaN(max_num_sv,N);
byLStemp  = NaN(max_num_sv,N);
H_posLStemp = NaN(max_num_sv,3,N);
nsvLS = NaN(1,N); 
dnsvLS = NaN(1,N); 
npriorLS = NaN(1,N); 
dnpriorLS = NaN(1,N);

err_LTS = NaN(p.LTSn,N);        hor_err_LTS = NaN(p.LTSn,N);  
ned_err_LTS = NaN(p.LTSn,N);    GDOPLTS = NaN(p.LTSn,N);
resLTS = cell(p.MShbn,1);       resLTStemp = NaN(max_num_sv*p.LTSn,N);
nsvLTS = NaN(p.LTSn,N);         dnsvLTS = NaN(p.LTSn,N); 
npriorLTS = NaN(p.LTSn,N);      dnpriorLTS = NaN(p.LTSn,N);
byLTStemp = NaN(max_num_sv*p.LTSn,N);
H_posLTStemp = NaN(max_num_sv*p.LTSn,3,N);

err_LTS2 = NaN(p.LTSn,N);       hor_err_LTS2 = NaN(p.LTSn,N);  
ned_err_LTS2 = NaN(p.LTSn,N);   GDOPLTS2 = NaN(p.LTSn,N);
resLTS2 = cell(p.MShbn,1);      resLTS2temp = NaN(max_num_sv*p.LTSn,N);
nsvLTS2 = NaN(p.LTSn,N);        dnsvLTS2 = NaN(p.LTSn,N); 
npriorLTS2 = NaN(p.LTSn,N);     dnpriorLTS2 = NaN(p.LTSn,N);
byLTS2temp = NaN(max_num_sv*p.LTSn,N);
H_posLTS2temp = NaN(max_num_sv*p.LTSn,3,N);

err_RAPS = NaN(p.RAPSn,N);      hor_err_RAPS = NaN(p.RAPSn,N); 
ned_err_RAPS = NaN(p.RAPSn,N);  GDOPRAPS = NaN(p.RAPSn,N);
resRAPS = cell(p.RAPSn,1);      resRAPStemp = NaN(max_num_sv*p.RAPSn,N);
nsvRAPS = NaN(p.RAPSn,N);       dnsvRAPS = NaN(p.RAPSn,N); 
npriorRAPS = NaN(p.RAPSn,N);    dnpriorRAPS= NaN(p.RAPSn,N);
byRAPStemp = NaN(max_num_sv*p.RAPSn,N);
H_posRAPStemp = NaN(max_num_sv*p.RAPSn,3,N);

err_TD = NaN(p.TDn,N);          hor_err_TD = NaN(p.TDn,N); 
ned_err_TD = NaN(p.TDn,N);      GDOPTD = NaN(p.TDn,N);
resTD = cell(p.TDn,1);          resTDtemp = NaN(max_num_sv*p.TDn,N);
nsvTD = NaN(p.TDn,N);           dnsvTD = NaN(p.TDn,N); 
npriorTD = NaN(p.TDn,N);        dnpriorTD= NaN(p.TDn,N);
byTDtemp = NaN(max_num_sv*p.TDn,N);
H_posTDtemp = NaN(max_num_sv*p.TDn,3,N);

err_TD2 = NaN(p.TDn,N);         hor_err_TD2 = NaN(p.TDn,N); 
ned_err_TD2 = NaN(p.TDn,N);     GDOPTD2 = NaN(p.TDn,N);
resTD2 = cell(p.TDn,1);         resTD2temp = NaN(max_num_sv*p.TDn,N);
nsvTD2 = NaN(p.TDn,N);          dnsvTD2 = NaN(p.TDn,N); 
npriorTD2 = NaN(p.TDn,N);       dnpriorTD2= NaN(p.TDn,N);
byTD2temp = NaN(max_num_sv*p.TDn,N);
H_posTD2temp = NaN(max_num_sv*p.TDn,3,N);

% For EKF -----------------------------
% state vector & covar matrix estimates at propagation steps
log.prior_state.LS       = nan(p.ns,p.numPropSteps*N); 
log.prior_statecovar.LS  = nan(p.ns,p.ns,p.numPropSteps*N);

% state vector & covar matrix estimates at measurement update steps
log.post_state.LS       = nan(p.ns,N+1); 
log.post_statecovar.LS  = nan(p.ns,p.ns,N+1);
log.delta_x.LS          = nan(p.ns,N+1);
log.Jydiag.LS               = nan(p.ns,N+1);

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

log.xhat = nan(13,N);

log.clkbNLS = nan(1,N+1);
log.msr_count = 0;
log.N = N;
got_new_xhat_posterior = true;

log.propTime = nan(1,N*p.numPropSteps);
log.postTime = nan(1,N+1); log.postTime(1) = 0;
for i = 1:N     
    p.idx = i;
    
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
    got_new_xhat_posterior = false;
    
    % Mark the sat prn that be computed
    gpslog_svprn_mark = []; glolog_svprn_mark = [];
    gallog_svprn_mark = []; bdslog_svprn_mark = [];
    % Record the prn of each system for whose satllite been used
    gpslog_prn_record = []; glolog_prn_record = [];
    gallog_prn_record = []; bdslog_prn_record = [];
    % Satellite position in ECEF frame
    gpslog_s_pos_ecef = []; glolog_s_pos_ecef = [];
    gallog_s_pos_ecef = []; bdslog_s_pos_ecef = [];
    % corrected pseudorange
    gpslog_corr_range = []; glolog_corr_range = [];
    gallog_corr_range = []; bdslog_corr_range = [];
    % The number of satellite be computed
    gpslog_num_sv = 0; glolog_num_sv = 0;
    gallog_num_sv = 0; bdslog_num_sv = 0;

    cpt = struct();
    if ~isempty(obs.GPS)&& p.freq==1 && p.enableGPS == 1
        % GPS satellite position computation, Single frenquency receiver mode    
        [gpslog_svprn_mark,gpslog_prn_record,gpslog_s_pos_ecef,gpslog_corr_range,gpslog_num_sv] =...
            satpost_corrpsedR_singlefreq_parfor(p,eph,obs,i,num_obs_gps,'GPS');
    end
    if ~isempty(obs.GLO)&& p.freq==1 && p.enableGLO == 1
        % GLO satellite position computation, Single frenquency receiver mode
        [glolog_svprn_mark,glolog_prn_record,glolog_s_pos_ecef,glolog_corr_range,glolog_num_sv] = ...
            satpost_corrpsedR_singlefreq_parfor(p,eph,obs,i,num_obs_glo,'GLO');
    end
    if ~isempty(obs.GAL)&& p.freq==1 && p.enableGAL == 1
        % GAL satellite position computation, Single frenquency receiver mode
        [gallog_svprn_mark,gallog_prn_record,gallog_s_pos_ecef,gallog_corr_range,gallog_num_sv] = ...
            satpost_corrpsedR_singlefreq_parfor(p,eph,obs,i,num_obs_gal,'GAL');
    end
    if ~isempty(obs.BDS)&& p.freq==1 && p.enableBDS == 1
        % BDS satellite position computation, Single frenquency receiver mode
        [bdslog_svprn_mark,bdslog_prn_record,bdslog_s_pos_ecef,bdslog_corr_range,bdslog_num_sv] = ...
            satpost_corrpsedR_singlefreq_parfor(p,eph,obs,i,num_obs_bds,'BDS');
    end
    cpt.prn_record = [gpslog_prn_record;glolog_prn_record;gallog_prn_record;bdslog_prn_record];
    cpt.svprn_mark = [gpslog_svprn_mark;glolog_svprn_mark;gallog_svprn_mark;bdslog_svprn_mark];
    cpt.corr_range = [gpslog_corr_range;glolog_corr_range;gallog_corr_range;bdslog_corr_range];
    ind = find(cpt.corr_range==0);
    cpt.corr_range(ind) = [];
    cpt.s_pos_ecef = [gpslog_s_pos_ecef,glolog_s_pos_ecef,gallog_s_pos_ecef,bdslog_s_pos_ecef];
    cpt.s_pos_ecef(:,ind) = [];
    cpt.num_sv = [gpslog_num_sv,glolog_num_sv,gallog_num_sv,bdslog_num_sv];
    % elevation & azimuth
    cpt.elev = NaN(length(cpt.corr_range),1); cpt.az = NaN(length(cpt.corr_range),1);
    % trop delay and iono delay
    cpt.trop_delay = NaN(length(cpt.corr_range),1); cpt.iono_delay = NaN(length(cpt.corr_range),1);
    
    if sum(cpt.num_sv)>=p.min_sv
        
        index = round(p.Grdpos.t) == round(obs.tr_sow(i));
        grdpos = p.Grdpos.pos(:,index);
        p.grdpos = grdpos; % DELETE ME !! after error check/debug  in each 
        % outlier accomodation algorithm solver is complete !!!!!!
        if isempty(grdpos)
            continue;
        end
        
        % Check elevation
        re_pos = p.x0(1:3);
        cpt = elevaz_check_linear(p,cpt,re_pos);
        if sum(cpt.num_sv)>=p.min_sv
            if ~isempty(p.eph_b) && ~isempty(p.obs_b)
                [cpt,n] = diff_corr_compute_linear(p,cpt,obs.tr_sow(i));
                if ~isempty(n) && sum(cpt.num_sv)>=p.min_sv
                    cpt = cpt_clear(cpt); % Clear the data where don't have diff correction
                    cpt.corr_range = cpt.corr_range - cpt.diff_corr;
                    %---------% Get prior for reciver clock bias
                    [~,x0_4,~] = userpos(p,cpt);
                    log.clkbNLS(i+1) = x0_4;
                    %---------%                    
                    [re_pos,clk_b,msr_res,cpt.GDOP,cpt.nsv,cpt.dnsv,cpt.nprior,cpt.dnprior,postxhat,postcovar,delta_x,Jydiag,by,H_pos] = ...
                        userposlinear_PVA(p,cpt.s_pos_ecef,cpt.corr_range,cpt.num_sv,x0_4,grdpos,x_prior,P_prior);
                    log.msr_count = log.msr_count + 1; % number of measurement updates
                    got_new_xhat_posterior = true;
                    log.postTime(i+1) = i; % GPS time when measurement update process is executed.
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
                    %--------------------save_result_linear---------------%
                    pos_ecef{i}   = re_pos;
                    rover_clk{i}  = clk_b;
                    sv_num_GPS(i) = cpt.num_sv(1);
                    sv_num_GLO(i) = cpt.num_sv(2);
                    sv_num_GAL(i) = cpt.num_sv(3);
                    sv_num_BDS(i) = cpt.num_sv(4);
                    ind_mark = cpt.svprn_mark ~= 0;
                    
                    
                    [ned_err_LS(:,i),hor_err_LS(:,i),err_LS(:,i),GDOPLS(:,i),...
                    nsvLS(i),dnsvLS(i),npriorLS(i),dnpriorLS(i),resLStemp(:,i),byLStemp(:,i),H_posLStemp(:,:,i)] = ...
                    save_errNorm_res_GDOP(1,grdpos,re_pos.LS,cpt.GDOP.LS,...
                    cpt.nsv.LS,cpt.dnsv.LS,cpt.nprior.LS,cpt.dnprior.LS,...
                    max_num_sv,ind_mark,msr_res.LS,by.LS,H_pos.LS);

                                      
                    if p.eb_LTS == 1
                        [ned_err_LTS(:,i),hor_err_LTS(:,i),err_LTS(:,i),GDOPLTS(:,i),...
                        nsvLTS(:,i),dnsvLTS(:,i),npriorLTS(:,i),dnpriorLTS(:,i),...
                        resLTStemp(:,i),byLTStemp(:,i),H_posLTStemp(:,:,i)] = save_errNorm_res_GDOP(p.LTSn,...
                        grdpos,re_pos.LTS,cpt.GDOP.LTS,cpt.nsv.LTS,cpt.dnsv.LTS,...
                        cpt.nprior.LTS,cpt.dnprior.LTS,max_num_sv,ind_mark,msr_res.LTS,by.LTS,H_pos.LTS);                                                
                    end
                    
                    if p.eb_RAPS == 1 
                        [ned_err_RAPS(:,i),hor_err_RAPS(:,i),err_RAPS(:,i),GDOPRAPS(:,i),...
                        nsvRAPS(:,i),dnsvRAPS(:,i),npriorRAPS(:,i),dnpriorRAPS(:,i),...
                        resRAPStemp(:,i),byRAPStemp(:,i),H_posRAPStemp(:,:,i)] = save_errNorm_res_GDOP(p.RAPSn,...
                        grdpos,re_pos.RAPS,cpt.GDOP.RAPS,cpt.nsv.RAPS,cpt.dnsv.RAPS,...
                        cpt.nprior.RAPS,cpt.dnprior.RAPS,max_num_sv,ind_mark,msr_res.RAPS,[],[]);                   
                    end
                    
                    if p.eb_TD == 1
                        [ned_err_TD(:,i),hor_err_TD(:,i),err_TD(:,i),GDOPTD(:,i),...
                        nsvTD(:,i),dnsvTD(:,i),npriorTD(:,i),dnpriorTD(:,i),...
                        resTDtemp(:,i),byTDtemp(:,i),H_posTDtemp(:,:,i)] = save_errNorm_res_GDOP(p.TDn,...
                        grdpos,re_pos.TD,cpt.GDOP.TD,cpt.nsv.TD,cpt.dnsv.TD,...
                        cpt.nprior.TD,cpt.dnprior.TD,max_num_sv,ind_mark,msr_res.TD,[],[]);                                           
                    end
                    
                    %-----------------------------------------------------%                    
                end
            else
                warning('No differential source given')
            end
       end
    end    
end

% output
log.tr_prime  = obs.tr_prime; % utc time
log.gpst      = obs.tr_sow-obs.tr_sow(1);
log.pos_ecef  = pos_ecef;     % Estimated position in ECEF
log.rover_clk = rover_clk;    % Receiver clock bias (meter)

log.sv_num_GPS = sv_num_GPS; 
log.sv_num_GLO = sv_num_GLO; 
log.sv_num_GAL = sv_num_GAL; 
log.sv_num_BDS = sv_num_BDS;

log.err_LS  = err_LS; log.hor_err_LS = hor_err_LS;
log.ned_err_LS = ned_err_LS;
log.GDOPLS = GDOPLS;
log.nsvLS = nsvLS;
log.dnsvLS = dnsvLS;
log.npriorLS = npriorLS;
log.dnpriorLS = dnpriorLS;
log.resLS = {resLStemp};
log.byLS = {byLStemp};
log.H_posLS = {H_posLStemp};

if p.eb_RAPS == 1
    pointer1 = 1;
    for idx = 1:p.RAPSn
        resRAPS{idx} = resRAPStemp(pointer1:pointer1+max_num_sv-1,:);
        byRAPS{idx} = byRAPStemp(pointer1:pointer1+max_num_sv-1,:);
        H_posRAPS{idx} = H_posRAPStemp(pointer1:pointer1+max_num_sv-1,3,:);
        pointer1 = pointer1 + max_num_sv;
    end
    log.err_RAPS = err_RAPS; log.hor_err_RAPS = hor_err_RAPS;
    log.ned_err_RAPS = ned_err_RAPS;
    log.GDOPRAPS = GDOPRAPS;
    log.nsvRAPS  = nsvRAPS;
    log.dnsvRAPS = dnsvRAPS;
    log.npriorRAPS = npriorRAPS;
    log.dnpriorRAPS = dnpriorRAPS;
    log.resRAPS = resRAPS;
    log.byRAPS  = byRAPS;
    log.H_posRAPS = H_posRAPS;
    log.RAPSEps = p.RAPSEps; log.RAPSn = p.RAPSn;
end


if p.eb_LTS == 1
    pointer1 = 1;
    for idx = 1:p.LTSn
        resLTS{idx} = resLTStemp(pointer1:pointer1+max_num_sv-1,:);
        byLTS{idx}  = byLTStemp(pointer1:pointer1+max_num_sv-1,:);
        H_posLTS{idx} = H_posLTStemp(pointer1:pointer1+max_num_sv-1,3,:);
        pointer1 = pointer1 + max_num_sv;
    end
    log.err_LTS = err_LTS; log.hor_err_LTS = hor_err_LTS; 
    log.ned_err_LTS = ned_err_LTS;
    log.GDOPLTS = GDOPLTS;
    log.nsvLTS  = nsvLTS;
    log.dnsvLTS = dnsvLTS;
    log.npriorLTS = npriorLTS;
    log.dnpriorLTS = dnpriorLTS;
    log.resLTS = resLTS;
    log.byLTS = byLTS;
    log.H_posLTS = H_posLTS;
    log.LTSOption = p.LTSOption; log.LTSn = p.LTSn;
end

if p.eb_TD == 1
    pointer1 = 1;
    for idx = 1:p.TDn
        resTD{idx} = resTDtemp(pointer1:pointer1+max_num_sv-1,:);
        byTD{idx}  = byTDtemp(pointer1:pointer1+max_num_sv-1,:);
        H_posTD{idx} = H_posTDtemp(pointer1:pointer1+max_num_sv-1,3,:);
        pointer1 = pointer1 + max_num_sv;
    end
    log.err_TD = err_TD; log.hor_err_TD = hor_err_TD; 
    log.ned_err_TD = ned_err_TD;
    log.GDOPTD = GDOPTD;
    log.nsvTD = nsvTD;
    log.dnsvTD = dnsvTD;
    log.npriorTD = npriorTD;
    log.dnpriorTD = dnpriorTD;
    log.resTD = resTD;
    log.byTD  = byTD;
    log.H_posTD = H_posTD;
    log.TDLambda = p.TDLambda; log.TDn = p.TDn;
end

log.p = p;
end