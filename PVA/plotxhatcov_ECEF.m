function [fig1,fig2,fig3] = plotxhatcov_ECEF(output,solvername,option)
    if solvername == "td",    TDn_i   = option.TDn_i;   end     % TDn_i   - ith TD estimation data
    if solvername == "lts",   LTSn_i  = option.LTSn_i;  end     % LTSn_i  - ith LTS estimation data
    if solvername == "raps",  RAPSn_i = option.RAPSn_i; end     % RAPSn_i - ith RAPS estimation data
    
    p = output.p;
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
    
    prpx_cov = x_prior_cov(1,1,:);  prvx_cov = x_prior_cov(4,4,:);  prax_cov = x_prior_cov(7,7,:);
    prpy_cov = x_prior_cov(2,2,:);  prvy_cov = x_prior_cov(5,5,:);  pray_cov = x_prior_cov(8,8,:);
    prpz_cov = x_prior_cov(3,3,:);  prvz_cov = x_prior_cov(6,6,:);  praz_cov = x_prior_cov(9,9,:);
    pspx_cov = x_post_cov(1,1,:);   psvx_cov = x_post_cov(4,4,:);   psax_cov = x_post_cov(7,7,:);
    pspy_cov = x_post_cov(2,2,:);   psvy_cov = x_post_cov(5,5,:);   psay_cov = x_post_cov(8,8,:);
    pspz_cov = x_post_cov(3,3,:);   psvz_cov = x_post_cov(6,6,:);   psaz_cov = x_post_cov(9,9,:);
    
    
    
    N = numel(output.gpst); % GPS epochs
    propsteps = output.p.numPropSteps;
%     total_propagation_steps = N*propsteps;
%     propagation_steps = (1:1:total_propagation_steps).*p.T;
    
    pltTstart1 = 1;    % <-------------------- plot time start     
    pltTend1   = 3750; % N; % <--------------- plot time end
    pltTstart2 = (pltTstart1-1) * propsteps + 1;
    pltTend2   = (pltTend1-1) * propsteps;    
    priorTime  = output.propTime(pltTstart2:pltTend2); %propagation_steps(pltTstart2:pltTend2);
    postTime   = output.gpst(pltTstart1:pltTend1);
    
    xlbl_txt = strcat('Receiver time using GPS second');
    xlimits = [0 3750];
    ylimits = [-inf inf];
    
    % Marker Properties
    mrkr = '.'; % marker
    sz = 10;    % marker size
    pr_mrk_clr = [0.55 0.55 0.55];
    ps_mrk_clr = 'r'; % marker color
    
    legend_txt = {'propagation estimate','posterior estimate'};
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
    figtag1 = strcat("cov_pos_",solvername);
    fignum = getfignum(figtag1);
    fig1 = figure(fignum); clf; hold on
    subplot(311); hold on; grid on
    scatter(priorTime,prpx_cov(pltTstart2:pltTend2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz)
    scatter(postTime,pspx_cov(pltTstart1:pltTend1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz)
    ylabel('cov(p_x),  m^2')
    xlim(xlimits); ylim(ylimits)
    legend(legend_txt);
    
    subplot(312); hold on; grid on
    scatter(priorTime,prpy_cov(pltTstart2:pltTend2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz)
    scatter(postTime,pspy_cov(pltTstart1:pltTend1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz)
    ylabel('cov(p_y),  m^2')
    xlim(xlimits); ylim(ylimits)
    legend(legend_txt);
    
    subplot(313); hold on; grid on
    scatter(priorTime,prpz_cov(pltTstart2:pltTend2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz)
    scatter(postTime,pspz_cov(pltTstart1:pltTend1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz)
    ylabel('cov(p_z),  m^2')
    xlim(xlimits); ylim(ylimits)
    xlabel(xlbl_txt)
    legend(legend_txt);
    sgtitle('Position Covariance in ECEF frame')
    
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------    
    figtag2 = strcat("cov_vel_",solvername);
    fignum = getfignum(figtag2);
    fig2 = figure(fignum); clf; hold on
    subplot(311); hold on; grid on
    scatter(priorTime,prvx_cov(pltTstart2:pltTend2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz)
    scatter(postTime,psvx_cov(pltTstart1:pltTend1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz)
    ylabel('cov(v_x),  m^2s^{-2}')
    xlim(xlimits); ylim(ylimits)
    legend(legend_txt);
    
    subplot(312); hold on; grid on
    scatter(priorTime,prvy_cov(pltTstart2:pltTend2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz)
    scatter(postTime,psvy_cov(pltTstart1:pltTend1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz)
    ylabel('cov(v_y),  m^2s^{-2}')
    xlim(xlimits); ylim(ylimits)
    legend(legend_txt);
    
    subplot(313); hold on; grid on
    scatter(priorTime,prvz_cov(pltTstart2:pltTend2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz)
    scatter(postTime,psvz_cov(pltTstart1:pltTend1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz)
    xlabel(xlbl_txt)
    ylabel('cov(v_z),  m^2s^{-2}')
    xlim(xlimits); ylim(ylimits)
    legend(legend_txt);
    sgtitle('Velocity Covariance in ECEF frame')
    
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
    figtag3 = strcat("cov_acl_",solvername);
    fignum = getfignum(figtag3);
    fig3 = figure(fignum); clf; hold on
    xlbl_txt = strcat('Receiver time using GPS second');
    xlbl_txt_msr = strcat('Receiver time using GPS second');
    subplot(311); hold on; grid on
    scatter(priorTime,prax_cov(pltTstart2:pltTend2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz)
    scatter(postTime,psax_cov(pltTstart1:pltTend1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz)
    ylabel('cov(a_x),  m^2s^{-4}')
    xlim(xlimits); ylim(ylimits)
    legend(legend_txt);
    
    subplot(312); hold on; grid on
    scatter(priorTime,pray_cov(pltTstart2:pltTend2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz)
    scatter(postTime,psay_cov(pltTstart1:pltTend1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz)
    ylabel('cov(a_y),  m^2s^{-4}')
    xlim(xlimits); ylim(ylimits)
    legend(legend_txt);
    
    subplot(313); hold on; grid on
    scatter(priorTime,praz_cov(pltTstart2:pltTend2),'Marker',mrkr,'MarkerEdgeColor',pr_mrk_clr,'SizeData',sz)
    scatter(postTime,psaz_cov(pltTstart1:pltTend1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz)
    xlabel(xlbl_txt)
    ylabel('cov(a_z),  m^2s^{-4}')
    xlim(xlimits); ylim(ylimits)
    legend(legend_txt);
    sgtitle('Acceleration Covariance in ECEF frame')    
    
end