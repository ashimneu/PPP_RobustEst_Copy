function [fig1,fig2,fig3] = plotxhat_NED(output,solvername,option)
    p = output.p;
    if solvername == "td",    TDn_i   = option.TDn_i;   end     % TDn_i   - ith TD estimation data
    if solvername == "lts",   LTSn_i  = option.LTSn_i;  end     % LTSn_i  - ith LTS estimation data
    if solvername == "raps",  RAPSn_i = option.RAPSn_i; end     % RAPSn_i - ith RAPS estimation data
    
    
    wgs84  = wgs84Ellipsoid(); % meters
    Grdpos = output.p.Grdpos.pos;  % Ground Truth pos in ECEF
    grd0   = ecef2lla(Grdpos(:,1)'); % [degrees] true initial pos
    [trueN,trueE,trueD] = ecef2ned(Grdpos(1,:),Grdpos(2,:),Grdpos(3,:),grd0(1),grd0(2),grd0(3),wgs84);
    
    % Origin of NED frame in ECEF frame, i.e. True Position at time = 0sec
    NED_origin = Grdpos(:,1);
    
    % Rotation matrix from ECEF to NED frame
    Rotation = [0.256585082665972   0.496560105722668   0.829211768114673;
                0.888404726250657  -0.459061043893141  -0.000000000413696;
                0.380658819413749   0.736675654127671  -0.558934561125126];    
    
    
    if lower(solvername) == "ls"
            x_prior = output.prior_state.LS;
            x_post  = output.post_state.LS;
    
    elseif lower(solvername) == "lts"
            x_prior = output.prior_state.LTS{LTSn_i};
            x_post  = output.post_state.LTS{LTSn_i};
            
    elseif lower(solvername) == "td"
            x_prior = output.prior_state.TD{TDn_i};
            x_post  = output.post_state.TD{TDn_i};
            
    elseif lower(solvername) == "raps"
            x_prior = output.prior_state.RAPS{RAPSn_i};
            x_post  = output.post_state.RAPS{RAPSn_i};        
    end

    [pos_prior(1,:),pos_prior(2,:),pos_prior(3,:)] = ecef2ned(x_prior(1,:),x_prior(2,:),x_prior(3,:),grd0(1),grd0(2),grd0(3),wgs84);
    [pos_post(1,:), pos_post(2,:), pos_post(3,:) ] = ecef2ned(x_post(1,:),x_post(2,:),x_post(3,:),grd0(1),grd0(2),grd0(3),wgs84);
    vel_prior = Rotation*x_prior(4:6,:);
    vel_post  = Rotation*x_post(4:6,:);            
    acl_prior = Rotation*x_prior(7:9,:);
    acl_post  = Rotation*x_post(7:9,:);            
    
    prpx = pos_prior(1,:);  prvx = vel_prior(1,:);  prax = acl_prior(1,:);
    prpy = pos_prior(2,:);  prvy = vel_prior(2,:);  pray = acl_prior(2,:);
    prpz = pos_prior(3,:);  prvz = vel_prior(3,:);  praz = acl_prior(3,:);
    pspx = pos_post(1,:);   psvx = vel_post(1,:);   psax = acl_post(1,:);
    pspy = pos_post(2,:);   psvy = vel_post(2,:);   psay = acl_post(2,:);
    pspz = pos_post(3,:);   psvz = vel_post(3,:);   psaz = acl_post(3,:);
    
    prct  = x_prior(10,:);   prcb = x_prior(11,:);
    psct  = x_post(10,:);    pscb = x_post(11,:);
    clkNLS = output.clkbNLS;
    
    N = size(Grdpos,2); % GPS epochs
    propsteps = output.p.numPropSteps;
%     total_propagation_steps = N*propsteps;
%     propagation_steps = (1:1:total_propagation_steps).*p.T;

    pltTstart1 = 1; % <-------------------- plot time start     
    pltTend1   = N; % <--------------- plot time end
    pltTstart2 = (pltTstart1-1) * propsteps + 1;
    pltTend2   = (pltTend1-1) * propsteps;
    priorTime  = output.propTime(pltTstart2:pltTend2); % priorTime = propagation_steps(pltTstart2:pltTend2);
    postTime   = output.gpst(pltTstart1:pltTend1); % begins at 0

    xlimits = [0 inf]; %[0 3770];  
    ylimits = [-inf inf];
    
    % Marker Properties
    mrkr = '.'; % marker
    sz = 10; % marker size
    tr_mrk_clr = 'b';
    pr_mrk_clr = [0.55 0.55 0.55];
    ps_mrk_clr = 'r'; % marker color
    
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

    figtag1 = strcat("pos_",solvername);
    fignum = getfignum(figtag1);    
    fig1 = figure(fignum); clf; hold on
    xlbl_txt = strcat('Receiver time using GPS second');
    legend_txt = {'true','propagation estimate','posterior estimate'};
    
    subplot(311); hold on;  hold on; grid on
    scatter(postTime,trueN(pltTstart1:pltTend1),'Marker',mrkr,'MarkerEdgeColor',tr_mrk_clr,'SizeData',sz);
    scatter(priorTime,prpx(pltTstart2:pltTend2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz);
    scatter(postTime,pspx(pltTstart1:pltTend1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz);    
    ylabel('p_N, m')
    xlim(xlimits); ylim(ylimits)
    legend(legend_txt);
       
    subplot(312); hold on; grid on
    scatter(postTime,trueE(pltTstart1:pltTend1),'Marker',mrkr,'MarkerEdgeColor',tr_mrk_clr,'SizeData',sz);
    scatter(priorTime,prpy(pltTstart2:pltTend2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz);
    scatter(postTime,pspy(pltTstart1:pltTend1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz);
    ylabel('p_E, m')
    xlim(xlimits); ylim(ylimits)
    legend(legend_txt);
    
    subplot(313); hold on; grid on
    scatter(postTime,trueD(pltTstart1:pltTend1),'Marker',mrkr,'MarkerEdgeColor',tr_mrk_clr,'SizeData',sz);
    scatter(priorTime,prpz(pltTstart2:pltTend2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz);
    scatter(postTime,pspz(pltTstart1:pltTend1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz); 
    ylabel('p_D, m')
    xlim(xlimits); ylim(ylimits)
    xlabel(xlbl_txt)
    sgtitle('Position in NED frame')    
    legend(legend_txt);
    
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------    
    
    figtag2 = strcat("vel_",solvername);
    fignum = getfignum(figtag2);
    fig2 = figure(fignum); clf; hold on
    xlbl_txt = strcat('Receiver time using GPS second');
    legend_txt = {'propagation estimate','posterior estimate'};
    
    subplot(311); hold on;  hold on; grid on
    scatter(priorTime,prvx(pltTstart2:pltTend2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz);
    scatter(postTime,psvx(pltTstart1:pltTend1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz); 
    ylabel('v_N, ms^{-1}')
    xlim(xlimits); ylim(ylimits)
    legend(legend_txt);
    
    subplot(312); hold on; grid on
    scatter(priorTime,prvy(pltTstart2:pltTend2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz);
    scatter(postTime,psvy(pltTstart1:pltTend1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz); 
    ylabel('v_E, ms^{-1}')
    xlim(xlimits); ylim(ylimits)
    legend(legend_txt);
    
    subplot(313); hold on; grid on
    scatter(priorTime,prvz(pltTstart2:pltTend2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz);
    scatter(postTime,psvz(pltTstart1:pltTend1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz); 
    xlabel(xlbl_txt)
    ylabel('v_D, ms^{-1}')
    xlim(xlimits); ylim(ylimits)
    sgtitle('Velocity in NED frame')
    legend(legend_txt);
    
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
    
    figtag3 = strcat("acl_",solvername);
    fignum = getfignum(figtag3);
    fig3 = figure(fignum); clf; hold on
    xlbl_txt = strcat('Receiver time using GPS second');
    
    subplot(311); hold on; grid on
    scatter(priorTime,prax(pltTstart2:pltTend2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz);
    scatter(postTime,psax(pltTstart1:pltTend1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz); 
    ylabel('a_N, ms^{-2}')
    xlim(xlimits); ylim(ylimits)
    legend(legend_txt);
    
    subplot(312); hold on; grid on
    scatter(priorTime,pray(pltTstart2:pltTend2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz);
    scatter(postTime,psay(pltTstart1:pltTend1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz); 
    ylabel('a_E, ms^{-2}')
    xlim(xlimits); ylim(ylimits)
    legend(legend_txt);
    
    subplot(313); hold on; grid on
    scatter(priorTime,praz(pltTstart2:pltTend2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz);
    scatter(postTime,psaz(pltTstart1:pltTend1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz); 
    xlabel(xlbl_txt)
    ylabel('a_D, ms^{-2}')
    xlim(xlimits); ylim(ylimits)
    sgtitle('Acceleration in NED frame')
    legend(legend_txt);
        
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
    
    figtag4 = strcat("clkdrift_",solvername);
    fignum = getfignum(figtag4);
    fig3 = figure(fignum); clf; hold on
    xlbl_txt = strcat('Receiver time using GPS second');
    legend_txt = {'true','propagation estimate','posterior estimate'};
    
    subplot(211); hold on; grid on
    scatter(postTime,clkNLS(pltTstart1:pltTend1),'Marker',mrkr,'MarkerEdgeColor',tr_mrk_clr,'SizeData',sz); 
    scatter(priorTime,prct(pltTstart2:pltTend2),'Marker', mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz);
    scatter(postTime,psct(pltTstart1:pltTend1),'Marker',  mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz); 
    ylabel('clock bias, m')
    xlim(xlimits); ylim(ylimits)
    legend(legend_txt);
    
    subplot(212); hold on; grid on
    scatter(priorTime,prcb(pltTstart2:pltTend2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz);
    scatter(postTime,pscb(pltTstart1:pltTend1),'Marker',mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz);
    ylabel('drift rate')
    xlabel(xlbl_txt)
    xlim(xlimits); ylim(ylimits)
    legend('propagation estimate','posterior estimate');
    
    sgtitle('Clock Bias & Drift Rate Estimates')
end