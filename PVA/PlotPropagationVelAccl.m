function [fig1,fig2,fig3] = PlotPropagationVelAccl(output,solvername,miscplot)
    
    if isempty(miscplot), miscplot = false; end
    clearfig = false; % if isempty(clearfig), clearfig = false; end
    fig2 = []; fig3 = [];
    p = output.p;
    
    if lower(solvername) == "ls"
        sgtitle_txt = strcat('KF (ECEF frame)');
    elseif lower(solvername) == "lts"
        sgtitle_txt = strcat('LTS (ECEF frame)');
    elseif lower(solvername) == "td"
        sgtitle_txt = strcat('TD (ECEF frame)');
    end
    
    N = numel(output.gpst); % GPS epochs
    propsteps = output.p.numPropSteps;
    total_propagation_steps = N*propsteps;
    propagation_steps = (1:1:total_propagation_steps).*p.T;
    
    pltTstart1 = 1;    % <-------------------- plot time start     
    pltTend1   = 3750; % N; % <--------------- plot time end
    pltTstart2 = (pltTstart1-1) * propsteps + 1;
    pltTend2   = (pltTend1-1) * propsteps;    
    priorTime = propagation_steps(pltTstart2:pltTend2);
    postTime  = output.gpst(pltTstart1:pltTend1);
    
    prior_vel = vecnorm(output.prior_state(4:6,:)); prior_hor_vel = vecnorm(output.prior_state(4:5,:));
    post_vel  = vecnorm(output.post_state(4:6,:));  post_hor_vel  = vecnorm(output.post_state(4:5,:));
    prior_acl = vecnorm(output.prior_state(7:9,:)); prior_hor_acl = vecnorm(output.prior_state(7:8,:));
    post_acl  = vecnorm(output.post_state(7:9,:));  post_hor_acl  = vecnorm(output.post_state(7:8,:));

    title_prop = 'Propagation';
    title_msr  = 'Msr Update';
    xlbl_txt_prop = strcat('Receiver time using GPS second');
    xlbl_txt_msr = strcat('Receiver time using GPS second');
    xlimits = [0 3750];
    ylimits = [-inf inf];
    sz = 10; % marker size
    Marker = '.';
    Marker_clr = [0.55 0.55 0.55];

    fig1 = figure; hold on
    if clearfig, clf; end 
    subplot(221)
    scatter(priorTime,prior_vel(pltTstart2:pltTend2),sz,Marker)
    title(title_prop)
    xlabel(xlbl_txt_prop)
    ylabel('3D velocity, unit: ms^{-1}');grid on
    xlim(xlimits); ylim(ylimits)
    
    subplot(222)
    scatter(postTime,post_vel(pltTstart1:pltTend1),sz,Marker)
    title(title_msr)
    xlabel(xlbl_txt_msr)
    ylabel('3D velocity, unit: ms^{-1}');grid on
    xlim(xlimits); ylim(ylimits)
    
    subplot(223)
    scatter(priorTime,prior_acl(pltTstart2:pltTend2),sz,Marker)
    title(title_prop)
    xlabel(xlbl_txt_prop)
    ylabel('3D acceleration, unit: ms^{-2}');grid on
    xlim(xlimits); ylim(ylimits)
    
    subplot(224)
    scatter(postTime,post_acl(pltTstart1:pltTend1),sz,Marker)
    title(title_msr)
    xlabel(xlbl_txt_msr)
    ylabel('3D acceleration, unit: ms^{-2}');grid on
    xlim(xlimits); ylim(ylimits)
    sgtitle(sgtitle_txt)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
    if miscplot
        fig2 = figure; hold on
        if clearfig, clf; end 
        xlbl_txt_prop = strcat('Propagation time step');
        subplot(221)
        scatter(priorTime,prior_hor_vel(pltTstart2:pltTend2),sz,Marker)
        title(title_prop)
        xlabel(xlbl_txt_prop)
        ylabel('hor. velocity, unit: ms^{-1}'); grid on
        xlim(xlimits); ylim(ylimits)
        
        subplot(222)
        scatter(priorTime,prior_vel(pltTstart2:pltTend2),sz,Marker)
        title(title_prop)
        xlabel(xlbl_txt_prop)
        ylabel('3D velocity, unit: ms^{-1}');grid on
        xlim(xlimits); ylim(ylimits)

        subplot(223)
        scatter(postTime,post_hor_vel(pltTstart1:pltTend1),sz,Marker)
        title(title_msr)
        xlabel(xlbl_txt_msr)
        ylabel('hor. velocity, unit: ms^{-1}'); grid on
        xlim(xlimits); ylim(ylimits)

        subplot(224)
        scatter(postTime,post_vel(pltTstart1:pltTend1),sz,Marker)
        title(title_msr)
        xlabel(xlbl_txt_msr)
        ylabel('3D velocity, unit: ms^{-1}');grid on
        xlim(xlimits); ylim(ylimits)
        sgtitle(sgtitle_txt)
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        fig3 = figure; hold on
        if clearfig, clf; end %#ok<*UNRCH>
        subplot(221)
        scatter(priorTime,prior_hor_acl(pltTstart2:pltTend2),sz,Marker)
        title(title_prop)
        xlabel(xlbl_txt_prop)
        ylabel('hor. acceleration, unit: ms^{-2}');grid on
        xlim(xlimits); ylim(ylimits)

        subplot(222)
        scatter(priorTime,prior_acl(pltTstart2:pltTend2),sz,Marker)
        title(title_prop)
        xlabel(xlbl_txt_prop)
        ylabel('3D acceleration, unit: ms^{-2}');grid on
        xlim(xlimits); ylim(ylimits)

        subplot(223)
        scatter(postTime,post_hor_acl(pltTstart1:pltTend1),sz,Marker)
        title(title_msr)
        xlabel(xlbl_txt_msr)
        ylabel('hor. acceleration, unit: ms^{-2}');grid on
        xlim(xlimits); ylim(ylimits)

        subplot(224)
        scatter(postTime,post_acl(pltTstart1:pltTend1),sz,Marker)
        title(title_msr)
        xlabel(xlbl_txt_msr)
        ylabel('3D acceleration, unit: ms^{-2}');grid on
        xlim(xlimits); ylim(ylimits)
        sgtitle(sgtitle_txt)
    end
end

