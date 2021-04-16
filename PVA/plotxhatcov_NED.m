function [fig1,fig2,fig3] = plotxhatcov_NED(output,solvername,option)
    p = output.p;
    if solvername == "td",    TDn_i   = option.TDn_i;   end     % TDn_i   - ith TD estimation data
    if solvername == "lts",   LTSn_i  = option.LTSn_i;  end     % LTSn_i  - ith LTS estimation data
    if solvername == "raps",  RAPSn_i = option.RAPSn_i; end     % RAPSn_i - ith RAPS estimation data
    
    
    wgs84  = wgs84Ellipsoid(); % meters
    Grdpos = output.p.Grdpos.pos;  % Ground Truth pos in ECEF
    p1 = Grdpos(:,1);
    p1_lla  = ecef2lla(Grdpos(:,1)'); % [degrees] true initial pos
    R = zeros(3,3); % Rotation matrix from ECEF to NED frame
    [R(1,1),R(2,1),R(3,1)] = ecef2ned(p1(1)+1,p1(2),p1(3),p1_lla(1),p1_lla(2),p1_lla(3),wgs84);
    [R(2,1),R(2,2),R(3,2)] = ecef2ned(p1(1),p1(2)+1,p1(3),p1_lla(1),p1_lla(2),p1_lla(3),wgs84);
    [R(1,3),R(2,3),R(3,3)] = ecef2ned(p1(1),p1(2),p1(3)+1,p1_lla(1),p1_lla(2),p1_lla(3),wgs84);
    
%     % Rotation matrix from ECEF to NED frame
%     Rotation = [0.256585082665972   0.496560105722668   0.829211768114673;
%                 0.888404726250657  -0.459061043893141  -0.000000000413696;
%                 0.380658819413749   0.736675654127671  -0.558934561125126];  
            
    if lower(solvername) == "ls"
        x_prior_cov = output.prior_statecovar.LS;
        x_post_cov  = output.post_statecovar.LS;
             
    elseif lower(solvername) == "lts"
        x_prior_cov = output.prior_statecovar.LTS{LTSn_i};
        x_post_cov  = output.post_statecovar.LTS{LTSn_i};
        
    elseif lower(solvername) == "td"
        x_prior_cov = output.prior_statecovar.TD{TDn_i};
        x_post_cov  = output.post_statecovar.TD{TDn_i};
        
    elseif lower(solvername) == "raps"
        x_prior_cov = output.prior_statecovar.RAPS{RAPSn_i};
        x_post_cov  = output.post_statecovar.RAPS{RAPSn_i};
          
    end
    
    % Convert from ECEF to NED frame
    pos_prior_cov = pagemtimes(pagemtimes(R,x_prior_cov(1:3,1:3,:)),R');
    pos_post_cov  = pagemtimes(pagemtimes(R,x_post_cov(1:3,1:3,:)),R');
    
    vel_prior_cov = pagemtimes(pagemtimes(R,x_prior_cov(4:6,4:6,:)),R');
    vel_post_cov  = pagemtimes(pagemtimes(R,x_post_cov(4:6,4:6,:)),R');
    
    acl_prior_cov = pagemtimes(pagemtimes(R,x_prior_cov(7:9,7:9,:)),R');
    acl_post_cov  = pagemtimes(pagemtimes(R,x_post_cov(7:9,7:9,:)),R');
    
    prpx_cov = pos_prior_cov(1,1,:);  prvx_cov = vel_prior_cov(1,1,:);  prax_cov = acl_prior_cov(1,1,:);
    prpy_cov = pos_prior_cov(2,2,:);  prvy_cov = vel_prior_cov(2,2,:);  pray_cov = acl_prior_cov(2,2,:);
    prpz_cov = pos_prior_cov(3,3,:);  prvz_cov = vel_prior_cov(3,3,:);  praz_cov = acl_prior_cov(3,3,:);
    pspx_cov = pos_post_cov(1,1,:);   psvx_cov = vel_post_cov(1,1,:);   psax_cov = acl_post_cov(1,1,:);
    pspy_cov = pos_post_cov(2,2,:);   psvy_cov = vel_post_cov(2,2,:);   psay_cov = acl_post_cov(2,2,:);
    pspz_cov = pos_post_cov(3,3,:);   psvz_cov = vel_post_cov(3,3,:);   psaz_cov = acl_post_cov(3,3,:);
    
    prct_cov = x_prior_cov(10,10,:); prcb_cov = x_prior_cov(11,11,:); 
    psct_cov = x_post_cov(10,10,:);  pscb_cov = x_post_cov(11,11,:);
    
    N = numel(output.gpst); % GPS epochs
    propsteps = output.p.numPropSteps;
%     total_propagation_steps = N*propsteps;
%     propagation_steps = (1:1:total_propagation_steps).*p.T;
    
    pltTstart1 = 1;    % <--------------- plot time start     
    pltTend1   = 3780; % <--------------- plot time end
    pltTstart2 = (pltTstart1-1) * propsteps + 1;
    pltTend2   = (pltTend1-1) * propsteps;    
    priorTime  = output.propTime(pltTstart2:pltTend2); % priorTime  = propagation_steps(pltTstart2:pltTend2);
    postTime   = output.gpst(pltTstart1:pltTend1);
    
    xlbl_txt = strcat('Receiver time using GPS second');
    xlimits = [0 3780];
    ylimits = [-inf inf];
    
    % Marker Properties
    mrkr = '.'; % marker
    sz = 10;    % marker size
    pr_mrk_clr = [0.55 0.55 0.55];
    ps_mrk_clr = 'r'; % marker color
    
    legend_txt = {'propagation estimate','posterior estimate'};
    
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
    figtag1 = strcat("cov_pos_",solvername); fignum = getfignum(figtag1);
    fig1 = figure(fignum); clf; hold on
    
    subplot(311); hold on; grid on
    scatter(priorTime,prpx_cov(pltTstart2:pltTend2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz)
    scatter(postTime,pspx_cov(pltTstart1:pltTend1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz)
    ylabel('cov(p_N), m^2'); 
    xlim(xlimits); ylim(ylimits)
    legend(legend_txt);
    
    subplot(312); hold on; grid on
    scatter(priorTime,prpy_cov(pltTstart2:pltTend2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz)
    scatter(postTime,pspy_cov(pltTstart1:pltTend1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz)
    ylabel('cov(p_E), m^2'); 
    xlim(xlimits); ylim(ylimits)
    legend(legend_txt);
    
    subplot(313); hold on; grid on
    scatter(priorTime,prpz_cov(pltTstart2:pltTend2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz)
    scatter(postTime,pspz_cov(pltTstart1:pltTend1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz)
    xlabel(xlbl_txt)
    ylabel('cov(p_D), m^2'); 
    xlim(xlimits); ylim(ylimits)
    legend(legend_txt);
    sgtitle('Position Covariance in NED frame')
        
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------    
    figtag2 = strcat("cov_vel_",solvername); fignum = getfignum(figtag2);
    fig2 = figure(fignum); clf; hold on
    
    subplot(311); hold on; grid on
    scatter(priorTime,prvx_cov(pltTstart2:pltTend2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz)
    scatter(postTime,psvx_cov(pltTstart1:pltTend1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz)
    ylabel('cov(v_N), m^2s^{-2}')
    xlim(xlimits); ylim(ylimits)
    legend(legend_txt);
    
    subplot(312); hold on; grid on
    scatter(priorTime,prvy_cov(pltTstart2:pltTend2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz)
    scatter(postTime,psvy_cov(pltTstart1:pltTend1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz)
    ylabel('cov(v_E), m^2s^{-2}')
    xlim(xlimits); ylim(ylimits)
    legend(legend_txt);
    
    subplot(313); hold on; grid on
    scatter(priorTime,prvz_cov(pltTstart2:pltTend2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz)
    scatter(postTime,psvz_cov(pltTstart1:pltTend1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz)
    xlabel(xlbl_txt)
    ylabel('cov(v_D), m^2s^{-2}')
    xlim(xlimits); ylim(ylimits)
    legend(legend_txt);
    sgtitle('Velocity Covariance in NED frame')
    
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
    figtag3 = strcat("cov_acl_",solvername); fignum = getfignum(figtag3);
    fig3 = figure(fignum); clf; hold on    

    subplot(311); hold on; grid on
    scatter(priorTime,prax_cov(pltTstart2:pltTend2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz)
    scatter(postTime,psax_cov(pltTstart1:pltTend1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz)
    ylabel('cov(a_N),  m^2s^{-4}')
    xlim(xlimits); ylim(ylimits)
    legend(legend_txt);
    
    subplot(312); hold on; grid on
    scatter(priorTime,pray_cov(pltTstart2:pltTend2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz)
    scatter(postTime,psay_cov(pltTstart1:pltTend1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz)
    ylabel('cov(a_E),  m^2s^{-4}')
    xlim(xlimits); ylim(ylimits)
    legend(legend_txt);
    
    subplot(313); hold on; grid on
    scatter(priorTime,praz_cov(pltTstart2:pltTend2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz)
    scatter(postTime,psaz_cov(pltTstart1:pltTend1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz)
    xlabel(xlbl_txt)
    ylabel('cov(a_D),  m^2s^{-4}')
    xlim(xlimits); ylim(ylimits)
    legend(legend_txt);
    sgtitle('Acceleration Covariance in NED frame')
    
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
    
    figtag4 = strcat("clkdrift_",solvername); fignum = getfignum(figtag4);
    fig4    = figure(fignum); clf; hold on
    xlbl_txt = strcat('Receiver time using GPS second');
    legend_txt = {'propagation estimate','posterior estimate'};
    
    subplot(211); hold on; grid on
    scatter(priorTime,prct_cov(pltTstart2:pltTend2),'Marker', mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz);
    scatter(postTime,psct_cov(pltTstart1:pltTend1),'Marker',  mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz); 
    ylabel('clock bias, m')
    xlim(xlimits); ylim(ylimits)
    legend(legend_txt);
    
    subplot(212); hold on; grid on
    scatter(priorTime,prcb_cov(pltTstart2:pltTend2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz);
    scatter(postTime,pscb_cov(pltTstart1:pltTend1),'Marker',mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz);
    ylabel('drift rate')
    xlabel(xlbl_txt)
    xlim(xlimits); ylim(ylimits)
    legend('propagation estimate','posterior estimate');
    
    sgtitle('Clock Bias & Drift Rate Covariances estimates')
    
end