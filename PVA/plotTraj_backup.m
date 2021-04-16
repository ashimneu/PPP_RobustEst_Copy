function option = plotTraj(output,solvername,option)
    p = output.p; ax = [];
    if solvername == "td",    TDn_i   = option.TDn_i;   end     % TDn_i   - ith TD estimation data
    if solvername == "lts",   LTSn_i  = option.LTSn_i;  end     % LTSn_i  - ith LTS estimation data
    if solvername == "raps",  RAPSn_i = option.RAPSn_i; end     % RAPSn_i - ith RAPS estimation data
    
    % Frame Conversion
    wgs84  = wgs84Ellipsoid(); % meters
    Grdpos = output.p.Grdpos.pos;  % Ground Truth pos in ECEF
    grd0   = ecef2lla(Grdpos(:,1)'); % [degrees] true initial pos
    [trueN,trueE,trueD] = ecef2ned(Grdpos(1,:),Grdpos(2,:),Grdpos(3,:),grd0(1),grd0(2),grd0(3),wgs84);
    
%     % Origin of NED frame in ECEF frame, i.e. True Position at time = 0sec
%     NED_origin = Grdpos(:,1);
%     p1_lla  = ecef2lla(Grdpos(:,1)'); % [degrees] true initial pos
%     T = zeros(3,3); % Rotation matrix from ECEF to NED frame
%     [T(1,1),T(2,1),T(3,1)] = ecef2ned(NED_origin(1)+1,NED_origin(2),NED_origin(3),p1_lla(1),p1_lla(2),p1_lla(3),wgs84);
%     [T(2,1),T(2,2),T(3,2)] = ecef2ned(NED_origin(1),NED_origin(2)+1,NED_origin(3),p1_lla(1),p1_lla(2),p1_lla(3),wgs84);
%     [T(1,3),T(2,3),T(3,3)] = ecef2ned(NED_origin(1),NED_origin(2),NED_origin(3)+1,p1_lla(1),p1_lla(2),p1_lla(3),wgs84);
    
    % Rotation matrix from ECEF to NED frame
    R = [0.256585082665972   0.496560105722668   0.829211768114673;
        0.888404726250657  -0.459061043893141  -0.000000000413696;
        0.380658819413749   0.736675654127671  -0.558934561125126];
    
    % Plot Time Window
    N = size(Grdpos,2); % GPS epochs
    propsteps = output.p.numPropSteps;

    pltTstart1 = 1;     % <--------------- plot time start     
    pltTend1   = N;  % <--------------- plot time end
    pltTstart2 = (pltTstart1-1) * propsteps + 1;
    pltTend2   = (pltTend1-1) * propsteps;
    TW1 = pltTstart1:pltTend1; % time window 1
    TW2 = pltTstart2:pltTend2; % time window 2
    priorTime = output.propTime(TW2); % priorTime = propagation_steps(pltTstart2:pltTend2);
    postTime  = output.gpst(TW1); % begins at 0
      
    
    if lower(solvername) == "kf"
        solvername = "kf";
        x_prior = output.prior_state.LS;
        x_post  = output.post_state.LS;
        x_prior_cov = output.prior_statecovar.LS;
        x_post_cov  = output.post_statecovar.LS;
        Jydiag = output.Jydiag.LS;
        msr_res = output.resLS{1};
        H_pos = output.H_posLS{1};
    
    elseif lower(solvername) == "lts"
        x_prior = output.prior_state.LTS{LTSn_i};
        x_post  = output.post_state.LTS{LTSn_i};
        x_prior_cov = output.prior_statecovar.LTS{LTSn_i};
        x_post_cov  = output.post_statecovar.LTS{LTSn_i};
        Jydiag = output.Jydiag.LS{LTSn_i};
        msr_res = output.resLTS(LTSn_i);
        H_pos = output.H_posLTS{LTSn_i};
            
    elseif lower(solvername) == "td"
        x_prior = output.prior_state.TD{TDn_i};
        x_post  = output.post_state.TD{TDn_i};
        x_prior_cov = output.prior_statecovar.TD{TDn_i};
        x_post_cov  = output.post_statecovar.TD{TDn_i};
        Jydiag = output.Jydiag.LS{TDn_i};
        msr_res = output.resLTS(TDn_i);
        H_pos = output.H_posTD{TDn_i};
            
    elseif lower(solvername) == "raps"
        x_prior = output.prior_state.RAPS{RAPSn_i};
        x_post  = output.post_state.RAPS{RAPSn_i};
        x_prior_cov = output.prior_statecovar.RAPS{RAPSn_i};
        x_post_cov  = output.post_statecovar.RAPS{RAPSn_i};
        Jydiag = output.Jydiag.LS{RAPSn_i};
        msr_res = output.resLTS(RAPSn_i);
        H_pos = output.H_posRAPS{RAPSn_i};
    end
    
    % Convert state estimates from ECEF to NED frame
    [pos_prior(1,:),pos_prior(2,:),pos_prior(3,:)] = ecef2ned(x_prior(1,:),x_prior(2,:),x_prior(3,:),grd0(1),grd0(2),grd0(3),wgs84);
    [pos_post(1,:), pos_post(2,:), pos_post(3,:) ] = ecef2ned(x_post(1,:),x_post(2,:),x_post(3,:),grd0(1),grd0(2),grd0(3),wgs84);
    vel_prior = R*x_prior(4:6,:);
    vel_post  = R*x_post(4:6,:);            
    acl_prior = R*x_prior(7:9,:);
    acl_post  = R*x_post(7:9,:);            
    
    prpx = pos_prior(1,:);  prvx = vel_prior(1,:);  prax = acl_prior(1,:);
    prpy = pos_prior(2,:);  prvy = vel_prior(2,:);  pray = acl_prior(2,:);
    prpz = pos_prior(3,:);  prvz = vel_prior(3,:);  praz = acl_prior(3,:);
    pspx = pos_post(1,:);   psvx = vel_post(1,:);   psax = acl_post(1,:);
    pspy = pos_post(2,:);   psvy = vel_post(2,:);   psay = acl_post(2,:);
    pspz = pos_post(3,:);   psvz = vel_post(3,:);   psaz = acl_post(3,:);
    
    cscale = 10^-5; % clock bias scale 
    cConst = output.clkbNLS(2); % constant value 
    clkNLS = (output.clkbNLS -cConst).*cscale; % NLS clk bias (true)
    
    prcb = (x_prior(10,:) - cConst).*cscale; % GPS Clock-Bias propagation
    pscb = (x_post(10,:) - cConst).*cscale;  % GPS Clock-Bias posterior
    prISBE = x_prior(11,:); % Inter-System Clock-Bias Galileo propagation
    psISBE = x_post(11,:);  % Inter-System Clock-Bias Galileo posterior
    prISBB = x_prior(12,:); % Inter-System Clock-Bias BeiDou propagation
    psISBB = x_post(12,:);  % Inter-System Clock-Bias BeiDou posterior
    prcd = x_prior(13,:);   % GPS Clock drift propagation      
    pscd = x_post(13,:);    % GPS Clock drift posterior   
    
    
    % Convert state covariances from ECEF to NED frame
    pos_prior_cov = pagemtimes(pagemtimes(R,x_prior_cov(1:3,1:3,:)),R');
    pos_post_cov  = pagemtimes(pagemtimes(R,x_post_cov(1:3,1:3,:)),R');    
    vel_prior_cov = pagemtimes(pagemtimes(R,x_prior_cov(4:6,4:6,:)),R');
    vel_post_cov  = pagemtimes(pagemtimes(R,x_post_cov(4:6,4:6,:)),R');    
    acl_prior_cov = pagemtimes(pagemtimes(R,x_prior_cov(7:9,7:9,:)),R');
    acl_post_cov  = pagemtimes(pagemtimes(R,x_post_cov(7:9,7:9,:)),R');
    
    prpx_cov = squeeze(pos_prior_cov(1,1,:));  prvx_cov   = squeeze(vel_prior_cov(1,1,:));  prax_cov   = squeeze(acl_prior_cov(1,1,:));
    prpy_cov = squeeze(pos_prior_cov(2,2,:));  prvy_cov   = squeeze(vel_prior_cov(2,2,:));  pray_cov   = squeeze(acl_prior_cov(2,2,:));
    prpz_cov = squeeze(pos_prior_cov(3,3,:));  prvz_cov   = squeeze(vel_prior_cov(3,3,:));  praz_cov   = squeeze(acl_prior_cov(3,3,:));
    pspx_cov = squeeze(pos_post_cov(1,1,:));   psvx_cov   = squeeze(vel_post_cov(1,1,:));   psax_cov   = squeeze(acl_post_cov(1,1,:));
    pspy_cov = squeeze(pos_post_cov(2,2,:));   psvy_cov   = squeeze(vel_post_cov(2,2,:));   psay_cov   = squeeze(acl_post_cov(2,2,:));
    pspz_cov = squeeze(pos_post_cov(3,3,:));   psvz_cov   = squeeze(vel_post_cov(3,3,:));   psaz_cov   = squeeze(acl_post_cov(3,3,:));
    prcb_cov = squeeze(x_prior_cov(10,10,:));  prISBE_cov = squeeze(x_prior_cov(11,11,:));  prISBB_cov = squeeze(x_prior_cov(12,12,:));
    pscb_cov = squeeze(x_post_cov(10,10,:));   psISBE_cov = squeeze(x_post_cov(11,11,:));   psISBB_cov = squeeze(x_post_cov(12,12,:));
    prcd_cov = squeeze(x_prior_cov(13,13,:));
    pscd_cov = squeeze(x_post_cov(13,13,:));
    
    prpx_std = sqrt(prpx_cov)';  prvx_std   = sqrt(prvx_cov)';   prax_std   = sqrt(prax_cov)';
    prpy_std = sqrt(prpy_cov)';  prvy_std   = sqrt(prvy_cov)';   pray_std   = sqrt(pray_cov)';
    prpz_std = sqrt(prpz_cov)';  prvz_std   = sqrt(prvz_cov)';   praz_std   = sqrt(praz_cov)';
    pspx_std = sqrt(pspx_cov)';  psvx_std   = sqrt(psvx_cov)';   psax_std   = sqrt(psax_cov)';
    pspy_std = sqrt(pspy_cov)';  psvy_std   = sqrt(psvy_cov)';   psay_std   = sqrt(psay_cov)';
    pspz_std = sqrt(pspz_cov)';  psvz_std   = sqrt(psvz_cov)';   psaz_std   = sqrt(psaz_cov)';
    prcb_std = sqrt(prcb_cov)'.*cscale;  
    pscb_std = sqrt(pscb_cov)'.*cscale;
    prISBE_std = sqrt(prISBE_cov)'; 
    psISBE_std = sqrt(psISBE_cov)';
    prISBB_std = sqrt(prISBB_cov)';
    psISBB_std = sqrt(psISBB_cov)';
    prcd_std = sqrt(prcd_cov)';
    pscd_std = sqrt(pscd_cov)';
    
    % Compute backdifferenced velocity: vel = (pos2 -pos1)/ (time2 - time1)
    % Here, backdiff. vel. is being referred to as true rover velocity.
    trueTime = output.p.Grdpos.t;
    trpx = trueN; trpy = trueE; trpz = trueD; 
    dt = trueTime(2:end) - trueTime(1:end-1); 
    % Padding with NaN because number of entries is 1 less after backdifferencing 
    trvx = [nan, (trpx(2:end) - trpx(1:end-1))./dt'];
    trvy = [nan, (trpy(2:end) - trpy(1:end-1))./dt'];
    trvz = [nan, (trpz(2:end) - trpz(1:end-1))./dt'];
    trvTime = postTime;

    % measurement residuals
    res = msr_res(any(msr_res,2),:); % get rows that aren't all nan
    
%--------------------------------------------------------------------------
% Figure paramters
%--------------------------------------------------------------------------
    xlimits = [0 inf]; %[0 3770];  
    ylimits = [-inf inf];
    
    % Marker Properties
    mrkr = '.'; % marker
    sz = 20; % marker size
    tr_mrk_clr = 'r';
    pr_mrk_clr = [0.55 0.55 0.55];
    ps_mrk_clr = 'b'; % marker color
    prstd_mrk_clr = 'c'; psstd_mrk_clr = 'c';
    
%--------------------------------------------------------------------------
% Position
%--------------------------------------------------------------------------

    figtag1 = strcat("pos_",solvername);
    fignum1 = getfignum(figtag1);    
    fig1 = figure(fignum1); clf; hold on
    xlbl_txt = strcat('Receiver time using GPS second');
    legend_txt = {'true','propagation est.','posterior est.','std'};
    
    subplot(311); hold on;  hold on; grid on
    scatter(postTime,trueN(TW1),'Marker',mrkr,'MarkerEdgeColor',tr_mrk_clr,'SizeData',sz);
    scatter(priorTime,prpx(TW2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz);
    scatter(priorTime,prpx(TW2)+ prpx_std(TW2),'Marker',mrkr,'MarkerEdgeColor',prstd_mrk_clr,'SizeData',sz,'HandleVisibility','off')
    scatter(priorTime,prpx(TW2)- prpx_std(TW2),'Marker',mrkr,'MarkerEdgeColor',prstd_mrk_clr,'SizeData',sz,'HandleVisibility','off')
    scatter(postTime,pspx(TW1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz);
    scatter(postTime,pspx(TW1)+pspx_std(TW1),'Marker', mrkr,'MarkerEdgeColor',psstd_mrk_clr,'SizeData',sz,'HandleVisibility','off');
    scatter(postTime,pspx(TW1)-pspx_std(TW1),'Marker', mrkr,'MarkerEdgeColor',psstd_mrk_clr,'SizeData',sz);
    ylabel('p_N, m')
    xlim(xlimits); ylim(ylimits)
           
    subplot(312); hold on; grid on
    scatter(postTime,trueE(TW1),'Marker',mrkr,'MarkerEdgeColor',tr_mrk_clr,'SizeData',sz);    
    scatter(priorTime,prpy(TW2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz);
    scatter(priorTime,prpy(TW2)+ prpy_std(TW2),'Marker',mrkr,'MarkerEdgeColor',prstd_mrk_clr,'SizeData',sz,'HandleVisibility','off')
    scatter(priorTime,prpy(TW2)- prpy_std(TW2),'Marker',mrkr,'MarkerEdgeColor',prstd_mrk_clr,'SizeData',sz,'HandleVisibility','off')
    scatter(postTime,pspy(TW1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz);
    scatter(postTime,pspy(TW1)+pspy_std(TW1),'Marker', mrkr,'MarkerEdgeColor',psstd_mrk_clr,'SizeData',sz,'HandleVisibility','off');
    scatter(postTime,pspy(TW1)-pspy_std(TW1),'Marker', mrkr,'MarkerEdgeColor',psstd_mrk_clr,'SizeData',sz);
    ylabel('p_E, m')
    xlim(xlimits); ylim(ylimits)
    
    splt = subplot(313); hold on; grid on
    scatter(postTime,trueD(TW1),'Marker',mrkr,'MarkerEdgeColor',tr_mrk_clr,'SizeData',sz);    
    scatter(priorTime,prpz(TW2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz);
    scatter(priorTime,prpz(TW2)+ prpz_std(TW2),'Marker',mrkr,'MarkerEdgeColor',prstd_mrk_clr,'SizeData',sz,'HandleVisibility','off')
    scatter(priorTime,prpz(TW2)- prpz_std(TW2),'Marker',mrkr,'MarkerEdgeColor',prstd_mrk_clr,'SizeData',sz,'HandleVisibility','off')
    scatter(postTime,pspz(TW1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz);
    scatter(postTime,pspz(TW1)+pspz_std(TW1),'Marker', mrkr,'MarkerEdgeColor',psstd_mrk_clr,'SizeData',sz,'HandleVisibility','off');
    scatter(postTime,pspz(TW1)-pspz_std(TW1),'Marker', mrkr,'MarkerEdgeColor',psstd_mrk_clr,'SizeData',sz);
    ylabel('p_D, m')
    xlim(xlimits); ylim(ylimits)
    xlabel(xlbl_txt)
    ax = [ax; findall(gcf, 'type', 'axes')];
    title_txt = strcat('Position in NED frame (',upper(solvername),')');
    sgtitle(title_txt)
    legend(legend_txt); adjustlegend(splt);
    
%--------------------------------------------------------------------------
% Velocity
%--------------------------------------------------------------------------    
    
    figtag2 = strcat("vel_",solvername);
    fignum2 = getfignum(figtag2);
    fig2 = figure(fignum2); clf; hold on
    xlbl_txt = strcat('Receiver time using GPS second');
    legend_txt = {'back diff.','propagation est.','posterior est.','std'};
    
    subplot(311); hold on;  hold on; grid on
    scatter(trvTime,trvx(TW1),'Marker',mrkr,'MarkerEdgeColor',tr_mrk_clr,'SizeData',sz); 
    scatter(priorTime,prvx(TW2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz);
    scatter(priorTime,prvx(TW2)+prvx_std(TW2),'Marker',mrkr,'MarkerEdgeColor',prstd_mrk_clr,'SizeData',sz,'HandleVisibility','off');
    scatter(priorTime,prvx(TW2)-prvx_std(TW2),'Marker',mrkr,'MarkerEdgeColor',prstd_mrk_clr,'SizeData',sz,'HandleVisibility','off');    
    scatter(postTime,psvx(TW1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz);
    scatter(postTime,psvx(TW1)+psvx_std(TW1),'Marker', mrkr,'MarkerEdgeColor',psstd_mrk_clr,'SizeData',sz,'HandleVisibility','off');
    scatter(postTime,psvx(TW1)-psvx_std(TW1),'Marker', mrkr,'MarkerEdgeColor',psstd_mrk_clr,'SizeData',sz);
    ylabel('v_N, ms^{-1}')
    xlim(xlimits); ylim(ylimits)
        
    subplot(312); hold on; grid on
    scatter(trvTime,trvy(TW1),'Marker',mrkr,'MarkerEdgeColor',tr_mrk_clr,'SizeData',sz); 
    scatter(priorTime,prvy(TW2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz);
    scatter(priorTime,prvy(TW2)+prvy_std(TW2),'Marker',mrkr,'MarkerEdgeColor',prstd_mrk_clr,'SizeData',sz,'HandleVisibility','off');
    scatter(priorTime,prvy(TW2)-prvy_std(TW2),'Marker',mrkr,'MarkerEdgeColor',prstd_mrk_clr,'SizeData',sz,'HandleVisibility','off');    
    scatter(postTime,psvy(TW1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz);
    scatter(postTime,psvy(TW1)+psvy_std(TW1),'Marker', mrkr,'MarkerEdgeColor',psstd_mrk_clr,'SizeData',sz,'HandleVisibility','off');
    scatter(postTime,psvy(TW1)-psvy_std(TW1),'Marker', mrkr,'MarkerEdgeColor',psstd_mrk_clr,'SizeData',sz);
    ylabel('v_E, ms^{-1}')
    xlim(xlimits); ylim(ylimits)
    
    splt = subplot(313); hold on; grid on
    scatter(trvTime,trvz(TW1),'Marker',mrkr,'MarkerEdgeColor',tr_mrk_clr,'SizeData',sz); 
    scatter(priorTime,prvz(TW2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz);
    scatter(priorTime,prvz(TW2)+prvz_std(TW2),'Marker',mrkr,'MarkerEdgeColor',prstd_mrk_clr,'SizeData',sz,'HandleVisibility','off');
    scatter(priorTime,prvz(TW2)-prvz_std(TW2),'Marker',mrkr,'MarkerEdgeColor',prstd_mrk_clr,'SizeData',sz,'HandleVisibility','off');    
    scatter(postTime,psvz(TW1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz);
    scatter(postTime,psvz(TW1)+psvz_std(TW1),'Marker', mrkr,'MarkerEdgeColor',psstd_mrk_clr,'SizeData',sz,'HandleVisibility','off');
    scatter(postTime,psvz(TW1)-psvz_std(TW1),'Marker', mrkr,'MarkerEdgeColor',psstd_mrk_clr,'SizeData',sz);
    xlabel(xlbl_txt)
    ylabel('v_D, ms^{-1}')
    xlim(xlimits); ylim(ylimits)
    ax = [ax; findall(gcf, 'type', 'axes')];
    title_txt = strcat('Velocity in NED frame (',upper(solvername),')');
    sgtitle(title_txt)
    legend(legend_txt); adjustlegend(splt);
    
%--------------------------------------------------------------------------
% Acceleration
%--------------------------------------------------------------------------
    
    figtag3 = strcat("acl_",solvername);
    fignum3 = getfignum(figtag3);
    fig3 = figure(fignum3); clf; hold on
    xlbl_txt = strcat('Receiver time using GPS second');
    legend_txt = {'propagation est.','posterior est.','std'};
    
    subplot(311); hold on; grid on
    scatter(priorTime,prax(TW2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz);
    scatter(priorTime,prax(TW2)+prax_std(TW2),'Marker',mrkr,'MarkerEdgeColor',prstd_mrk_clr,'SizeData',sz,'HandleVisibility','off');
    scatter(priorTime,prax(TW2)-prax_std(TW2),'Marker',mrkr,'MarkerEdgeColor',prstd_mrk_clr,'SizeData',sz,'HandleVisibility','off');    
    scatter(postTime,psax(TW1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz);
    scatter(postTime,psax(TW1)+psax_std(TW1),'Marker', mrkr,'MarkerEdgeColor',psstd_mrk_clr,'SizeData',sz,'HandleVisibility','off');
    scatter(postTime,psax(TW1)-psax_std(TW1),'Marker', mrkr,'MarkerEdgeColor',psstd_mrk_clr,'SizeData',sz);
    ylabel('a_N, ms^{-2}')
    xlim(xlimits); ylim(ylimits)
        
    subplot(312); hold on; grid on
    scatter(priorTime,pray(TW2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz);
    scatter(priorTime,pray(TW2)+pray_std(TW2),'Marker',mrkr,'MarkerEdgeColor',prstd_mrk_clr,'SizeData',sz,'HandleVisibility','off');
    scatter(priorTime,pray(TW2)-pray_std(TW2),'Marker',mrkr,'MarkerEdgeColor',prstd_mrk_clr,'SizeData',sz,'HandleVisibility','off');    
    scatter(postTime,psay(TW1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz);
    scatter(postTime,psay(TW1)+psay_std(TW1),'Marker', mrkr,'MarkerEdgeColor',psstd_mrk_clr,'SizeData',sz,'HandleVisibility','off');
    scatter(postTime,psay(TW1)-psay_std(TW1),'Marker', mrkr,'MarkerEdgeColor',psstd_mrk_clr,'SizeData',sz);
    ylabel('a_E, ms^{-2}')
    xlim(xlimits); ylim(ylimits)
    
    splt = subplot(313); hold on; grid on
    scatter(priorTime,praz(TW2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz);
    scatter(priorTime,praz(TW2)+praz_std(TW2),'Marker',mrkr,'MarkerEdgeColor',prstd_mrk_clr,'SizeData',sz,'HandleVisibility','off');
    scatter(priorTime,praz(TW2)-praz_std(TW2),'Marker',mrkr,'MarkerEdgeColor',prstd_mrk_clr,'SizeData',sz,'HandleVisibility','off');  
    scatter(postTime,psaz(TW1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz);
    scatter(postTime,psaz(TW1)+psaz_std(TW1),'Marker', mrkr,'MarkerEdgeColor',psstd_mrk_clr,'SizeData',sz,'HandleVisibility','off');
    scatter(postTime,psaz(TW1)-psaz_std(TW1),'Marker', mrkr,'MarkerEdgeColor',psstd_mrk_clr,'SizeData',sz);
    xlabel(xlbl_txt)
    ylabel('a_D, ms^{-2}')
    xlim(xlimits); ylim(ylimits)
    ax = [ax; findall(gcf, 'type', 'axes')];
    title_txt = strcat('Acceleration in NED frame (',upper(solvername),')');
    sgtitle(title_txt)
    legend(legend_txt); adjustlegend(splt);
        
%--------------------------------------------------------------------------
% Clock Drift
%--------------------------------------------------------------------------
    figtag4 = strcat("clk_biasdrift_",solvername);
    fignum4 = getfignum(figtag4);
    fig4 = figure(fignum4); clf; hold on
    xlbl_txt = strcat('Receiver time using GPS second');
    legend_txt = {'true','propagation est.','posterior est.','std'};
    
    splt1 = subplot(211); hold on; grid on
    scatter(postTime,clkNLS(TW1),'Marker',mrkr,'MarkerEdgeColor',tr_mrk_clr,'SizeData',sz); 
    
    scatter(priorTime,prcb(TW2),'Marker', mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz);
    scatter(priorTime,prcb(TW2)+prcb_std(TW2),'Marker',mrkr,'MarkerEdgeColor',prstd_mrk_clr,'SizeData',sz,'HandleVisibility','off');
    scatter(priorTime,prcb(TW2)-prcb_std(TW2),'Marker',mrkr,'MarkerEdgeColor',prstd_mrk_clr,'SizeData',sz,'HandleVisibility','off');    
    scatter(postTime,pscb(TW1),'Marker',  mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz); 
    scatter(postTime,pscb(TW1)+pscb_std(TW1),'Marker',mrkr,'MarkerEdgeColor',psstd_mrk_clr,'SizeData',sz,'HandleVisibility','off')
    scatter(postTime,pscb(TW1)-pscb_std(TW1),'Marker',mrkr,'MarkerEdgeColor',psstd_mrk_clr,'SizeData',sz);
    ylabel('clock bias $\times 10^{5}$, m','Interpreter','latex')
    xlim(xlimits); ylim(ylimits)
    legend(legend_txt); adjustlegend(splt1,'bottomleft');
    
    splt2 = subplot(212); hold on; grid on
    scatter(priorTime,prcd(TW2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz);
    scatter(priorTime,prcd(TW2)+prcd_std(TW2),'Marker',mrkr,'MarkerEdgeColor',prstd_mrk_clr,'SizeData',sz,'HandleVisibility','off');
    scatter(priorTime,prcd(TW2)-prcd_std(TW2),'Marker',mrkr,'MarkerEdgeColor',prstd_mrk_clr,'SizeData',sz,'HandleVisibility','off');    
    scatter(postTime,pscd(TW1),'Marker',mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz);
    scatter(postTime,pscd(TW1)+pscd_std(TW1),'Marker',mrkr,'MarkerEdgeColor',psstd_mrk_clr,'SizeData',sz,'HandleVisibility','off')
    scatter(postTime,pscd(TW1)-pscd_std(TW1),'Marker',mrkr,'MarkerEdgeColor',psstd_mrk_clr,'SizeData',sz);
    ylabel('drift rate')
    xlabel(xlbl_txt)
    xlim(xlimits); ylim(ylimits)
    
%     legend(legend_txt(2:end)); adjustlegend(splt2);
    ax = [ax; findall(gcf, 'type', 'axes')];
    title_txt = strcat('Clock Bias & Drift (',upper(solvername),')');
    sgtitle(title_txt)
    
%--------------------------------------------------------------------------
% Inter-System Clock Bias
%--------------------------------------------------------------------------

    figtag5 = strcat("ISB_",solvername);
    fignum5 = getfignum(figtag5);
    fig5 = figure(fignum5); clf; hold on
    xlbl_txt = strcat('Receiver time using GPS second');
    legend_txt = {'propagation est.','posterior est.','std'};
    
    subplot(211); hold on; grid on    
    scatter(priorTime,prISBE(TW2),'Marker', mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz);
    scatter(priorTime,prISBE(TW2)+prISBE_std(TW2),'Marker',mrkr,'MarkerEdgeColor',prstd_mrk_clr,'SizeData',sz,'HandleVisibility','off');
    scatter(priorTime,prISBE(TW2)-prISBE_std(TW2),'Marker',mrkr,'MarkerEdgeColor',prstd_mrk_clr,'SizeData',sz,'HandleVisibility','off');    
    scatter(postTime,psISBE(TW1),'Marker',  mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz); 
    scatter(postTime,psISBE(TW1)+psISBE_std(TW1),'Marker',mrkr,'MarkerEdgeColor',psstd_mrk_clr,'SizeData',sz,'HandleVisibility','off')
    scatter(postTime,psISBE(TW1)-psISBE_std(TW1),'Marker',mrkr,'MarkerEdgeColor',psstd_mrk_clr,'SizeData',sz);
    ylabel('ISB Galileo, m')
    xlim(xlimits); ylim(ylimits)
    
    splt = subplot(212); hold on; grid on
    scatter(priorTime,prISBB(TW2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz);
    scatter(priorTime,prISBB(TW2)+prISBB_std(TW2),'Marker',mrkr,'MarkerEdgeColor',prstd_mrk_clr,'SizeData',sz,'HandleVisibility','off');
    scatter(priorTime,prISBB(TW2)-prISBB_std(TW2),'Marker',mrkr,'MarkerEdgeColor',prstd_mrk_clr,'SizeData',sz,'HandleVisibility','off');    
    scatter(postTime,psISBB(TW1),'Marker',mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz);
    scatter(postTime,psISBB(TW1)+psISBB_std(TW1),'Marker',mrkr,'MarkerEdgeColor',psstd_mrk_clr,'SizeData',sz,'HandleVisibility','off')
    scatter(postTime,psISBB(TW1)-psISBB_std(TW1),'Marker',mrkr,'MarkerEdgeColor',psstd_mrk_clr,'SizeData',sz);
    ylabel('ISB BeiDou, m')
    xlabel(xlbl_txt)
    xlim(xlimits); ylim(ylimits)    
    
    legend(legend_txt); adjustlegend(splt);    
    ax = [ax; findall(gcf, 'type', 'axes')];
    title_txt = strcat('Inter-System Clock Bias (ISB) (',upper(solvername),')');
    sgtitle(title_txt)
    
%--------------------------------------------------------------------------
% Measurement Residuals
%--------------------------------------------------------------------------

    figtag6 = strcat("res_",solvername);
    fignum6 = getfignum(figtag6);
    fig6 = figure(fignum6); clf; hold on; grid on
    xlbl_txt = strcat('Receiver time using GPS second');    
    
    for i=1:1:size(res,1)   
        scatter(postTime,res(i,TW1),'Marker',mrkr,'SizeData',sz,'HandleVisibility','off'); 
    end
    xlabel(xlbl_txt)
    ylabel('residual (m)')
    ax = [ax; findall(gcf, 'type', 'axes')];
    title_txt = strcat('Measurement Residual (',upper(solvername),')');
    sgtitle(title_txt)
        
%--------------------------------------------------------------------------
% Information
%--------------------------------------------------------------------------
% if option.frame doesn't exist, set the default.
    if ~isfield(option, 'frame')
        option.frame = "ned"; % default frame
    end
    
    if option.frame == "ecef"
        ifpx = Jydiag(1,:); 
        ifpy = Jydiag(2,:); 
        ifpz = Jydiag(3,:);
    elseif option.frame == "ned"
        % Given rotation matrix,Rot, this section transforms Measurement 
        % Information (i.e. J_y = H'*Cov^(-1)*H) which is originally 
        % in ECEF frame to its representation in NED frame 
        % ie. Rot'*H'*Cov^(-1)*H*Rot

        J_pos = nan(3,3,size(H_pos,3));
        for page = 1:size(H_pos,3)
            H_nan = H_pos(:,:,page);    % H contains rows of nan
            H  = H_nan(any(H_nan,2),:); % rows of nan are removed
            m  = size(H,1)/2;
            invYCov = eye(m).*output.p.sig_y^2; % measurement noise covariance
            J_pos(:,:,page) = R'*H(1:m,1:3)'*invYCov*H*R;
        end        
        ifpx   = squeeze(J_pos(1,1,:)); 
        ifpy   = squeeze(J_pos(2,2,:)); 
        ifpz   = squeeze(J_pos(3,3,:));
    end
    ifcb   = Jydiag(10,:);
    ifISBE = Jydiag(11,:); 
    ifISBB = Jydiag(12,:);  
%--------------------------------------------------------------------------    
    figtag7 = strcat("info_pos_",solvername);
    fignum7 = getfignum(figtag7);    
    fig7 = figure(fignum7); clf; hold on
    xlbl_txt = strcat('Receiver time using GPS second');
    legend_txt = {'info.'};  

    subplot(311); hold on; grid on
    scatter(postTime,ifpx(TW1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz);
    ylabel('J_{p_N}, m^{-2}')
    xlim(xlimits); ylim(ylimits)

    subplot(312); hold on; grid on 
    scatter(postTime,ifpy(TW1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz);
    ylabel('J_{p_E}, m^{-2})')
    xlim(xlimits); ylim(ylimits)

    splt = subplot(313); hold on; grid on    
    scatter(postTime,ifpz(TW1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz);
    xlabel(xlbl_txt)
    ylabel('J_{p_D}, m^{-2}')
    xlim(xlimits); ylim(ylimits)
    ax = [ax; findall(gcf, 'type', 'axes')];
    title_txt = strcat('Information from Measurements in NED frame (',upper(solvername),')');
    sgtitle(title_txt)
    legend(legend_txt); adjustlegend(splt);

%--------------------------------------------------------------------------
% Clock Information
%--------------------------------------------------------------------------
            
    figtag8 = strcat("info_clk_",solvername);
    fignum8 = getfignum(figtag8);
    fig8 = figure(fignum8); clf; hold on
    xlbl_txt = strcat('Receiver time using GPS second');
    legend_txt = {'information'};
    
    subplot(311); hold on; grid on
    scatter(postTime,ifcb(TW1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz);
    ylabel('J_{ctr}, ms^{-2}')
    xlim(xlimits); ylim(ylimits)
           
    subplot(312); hold on; grid on 
    scatter(postTime,ifISBE(TW1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz);
    ylabel('J_{ISB_E}, ms^{-2}')
    xlim(xlimits); ylim(ylimits)
    
    splt = subplot(313); hold on; grid on    
    scatter(postTime,ifISBB(TW1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz);
    xlabel(xlbl_txt)
    ylabel('J_{ISB_B}, ms^{-2}')
    xlim(xlimits); ylim(ylimits)
    ax = [ax; findall(gcf, 'type', 'axes')];
    title_txt = strcat('Clock Information (',upper(solvername),')');
    sgtitle(title_txt)
    legend(legend_txt); adjustlegend(splt); 

%--------------------------------------------------------------------------
% Position Error
%--------------------------------------------------------------------------
    
    figtag9 = strcat("pos_err_",solvername);
    fignum9 = getfignum(figtag9);    
    fig9 = figure(fignum9); clf; hold on
    xlbl_txt = strcat('Receiver time using GPS second');
    legend_txt = {'estimation error','std'};
    ned = 1; ecef = 0;
    if ned
        er_pspx = pspx(1:N) - trueN(1:N);
        er_pspy = pspy(1:N) - trueE(1:N);
        er_pspz = pspz(1:N) - trueD(1:N);
    else
        er_pspx = pspx(1:N) - trueN(1:N);
        er_pspy = pspy(1:N) - trueE(1:N);
        er_pspz = pspz(1:N) - trueD(1:N);
    end
    
    subplot(311); hold on;  hold on; grid on
    scatter(postTime,er_pspx(TW1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz,'HandleVisibility','off');
    scatter(postTime,pspx_std(TW1),'Marker', mrkr,'MarkerEdgeColor',psstd_mrk_clr,'SizeData',sz,'HandleVisibility','off');
    scatter(postTime,-pspx_std(TW1),'Marker', mrkr,'MarkerEdgeColor',psstd_mrk_clr,'SizeData',sz,'HandleVisibility','off');
    ylabel('p_N, m')
    xlim(xlimits); ylim(ylimits)
           
    subplot(312); hold on; grid on
    scatter(postTime,er_pspy(TW1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz,'HandleVisibility','off');
    scatter(postTime,pspy_std(TW1),'Marker', mrkr,'MarkerEdgeColor',psstd_mrk_clr,'SizeData',sz,'HandleVisibility','off');
    scatter(postTime,-pspy_std(TW1),'Marker', mrkr,'MarkerEdgeColor',psstd_mrk_clr,'SizeData',sz,'HandleVisibility','off');
    ylabel('p_E, m')
    xlim(xlimits); ylim(ylimits)
    
    splt = subplot(313); hold on; grid on
    scatter(postTime,er_pspz(TW1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz);
    scatter(postTime,pspz_std(TW1),'Marker', mrkr,'MarkerEdgeColor',psstd_mrk_clr,'SizeData',sz);
    scatter(postTime,-pspz_std(TW1),'Marker', mrkr,'MarkerEdgeColor',psstd_mrk_clr,'SizeData',sz,'HandleVisibility','off');
    ylabel('p_D, m')
    xlim(xlimits); ylim(ylimits)
    xlabel(xlbl_txt)
    ax = [ax; findall(gcf, 'type', 'axes')];
    title_txt = strcat('Position error in NED frame (',upper(solvername),')');
    sgtitle(title_txt)
    legend(legend_txt); adjustlegend(splt);
    
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

    % Link axis limits of all the figures
    if ~isempty(ax) 
        if isfield(option,'ax')
            option.ax = [option.ax; ax];
        else
            option.ax = ax;
        end
    end
    
%     if isfield(option,'axes2link')
%         if ~strcmp(option.axes2link,'') | ~isempty(option.axes2link)
%             linkaxes(ax,option.axes2link)
%         end
%     end
end