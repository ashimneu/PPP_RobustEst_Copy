function [fig1,fig4] = plotInfo(output,solvername,option)
    p = output.p; ax = [];
    
    % Convert frame from ECEF to NED
    wgs84  = wgs84Ellipsoid(); % meters
    Grdpos = output.p.Grdpos.pos;  % Ground Truth pos in ECEF
    grd0   = ecef2lla(Grdpos(:,1)'); % [degrees] true initial pos
    % true position in NED frame
    [trueN,trueE,trueD] = ecef2ned(Grdpos(1,:),Grdpos(2,:),Grdpos(3,:),grd0(1),grd0(2),grd0(3),wgs84);  
    
%     % Origin of NED frame in ECEF frame, i.e. True Position at time = 0sec
%     NED_origin = Grdpos(:,1);
%     p1_lla  = ecef2lla(Grdpos(:,1)'); % [degrees] true initial pos
%     T = zeros(3,3); % Rotation matrix from ECEF to NED frame
%     [T(1,1),T(2,1),T(3,1)] = ecef2ned(NED_origin(1)+1,NED_origin(2),NED_origin(3),p1_lla(1),p1_lla(2),p1_lla(3),wgs84);
%     [T(2,1),T(2,2),T(3,2)] = ecef2ned(NED_origin(1),NED_origin(2)+1,NED_origin(3),p1_lla(1),p1_lla(2),p1_lla(3),wgs84);
%     [T(1,3),T(2,3),T(3,3)] = ecef2ned(NED_origin(1),NED_origin(2),NED_origin(3)+1,p1_lla(1),p1_lla(2),p1_lla(3),wgs84);
    
    % Rotation matrix from ECEF to NED frame
    Rot = [0.256585082665972   0.496560105722668   0.829211768114673;
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
        Jydiag = output.Jydiag.LS;
        H_pos = output.H_posLS{1};
    
    elseif lower(solvername) == "lts"
        LTSn_i  = option.LTSn_i;
        Jydiag = output.Jydiag.LS{LTSn_i};
        H_pos = output.H_posLTS{LTSn_i};
            
    elseif lower(solvername) == "td"
        TDn_i   = option.TDn_i;
        Jydiag = output.Jydiag.LS{TDn_i};
        H_pos = output.H_posTD{TDn_i};
            
    elseif lower(solvername) == "raps"
        RAPSn_i = option.RAPSn_i;
        Jydiag = output.Jydiag.LS{RAPSn_i};
        H_pos = output.H_posRAPS{RAPSn_i};
    end
    
    % if option.frame doesn't exist, set the default.
    if ~isfield(option, 'frame')
        option.frame = "ned"; % default frame
    end
    
    if option.frame == "ecef"
        ifpx = Jydiag(1,:); 
        ifpy = Jydiag(2,:); 
        ifpz = Jydiag(3,:);
    elseif option.frame == "ned"
        % Given rotation matrix,Rot, this section transforms Measurement Information (i.e. J_y = H'*Cov^(-1)*H)
        % which is originally in ECEF frame to its representation in NED 
        % frame ie. Rot'*H'*Cov^(-1)*H*Rot

        J_pos = nan(3,3,size(H_pos,3));
        for page = 1:size(H_pos,3)
            H_nan = H_pos(:,:,page);    % H contains rows of nan
            H  = H_nan(any(H_nan,2),:); % rows of nan are removed
            m  = size(H,1);
            iR = eye(m).*output.p.sig_y^2; % noise covariance
            J_pos(:,:,page) = Rot'*H'*iR*H*Rot;
        end        
        ifpx   = squeeze(J_pos(1,1,:)); 
        ifpy   = squeeze(J_pos(2,2,:)); 
        ifpz   = squeeze(J_pos(3,3,:));
    end
    ifcb   = Jydiag(10,:);
    ifISBE = Jydiag(11,:); 
    ifISBB = Jydiag(12,:);
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
% Figures
%--------------------------------------------------------------------------
    if option.frame == "ecef"
                
        figtag1 = strcat("info_pos_ECEF_",solvername);
        fignum1 = getfignum(figtag1);    
        fig1 = figure(fignum1); clf; hold on
        xlbl_txt = strcat('Receiver time using GPS second');
        legend_txt = {'info.'};  

        subplot(311); hold on; grid on
        scatter(postTime,ifpx(TW1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz);
        ylabel('J_{p_x}, m^{-2}')
        xlim(xlimits); ylim(ylimits)

        subplot(312); hold on; grid on 
        scatter(postTime,ifpy(TW1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz);
        ylabel('J_{p_y}, m^{-2})')
        xlim(xlimits); ylim(ylimits)

        splt = subplot(313); hold on; grid on    
        scatter(postTime,ifpz(TW1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz);
        xlabel(xlbl_txt)
        ylabel('J_{p_z}, m^{-2}')
        xlim(xlimits); ylim(ylimits)
        ax = [ax; findall(gcf, 'type', 'axes')];
        title_txt = strcat('Position Information in ECEF frame, (',upper(solvername),')');
        sgtitle(title_txt)
        legend(legend_txt); adjustlegend(splt);

    elseif option.frame == "ned"
        
        figtag1 = strcat("info_pos_NED_",solvername);
        fignum1 = getfignum(figtag1);    
        fig1 = figure(fignum1); clf; hold on
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
        title_txt = strcat('Position Information in NED frame, (',upper(solvername),')');
        sgtitle(title_txt)
        legend(legend_txt); adjustlegend(splt);
    end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
            
    figtag4 = strcat("info_clock_",solvername);
    fignum4 = getfignum(figtag4);
    fig4 = figure(fignum4); clf; hold on
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
    title_txt = strcat('Clock Information, (',upper(solvername),')');
    sgtitle(title_txt)
    legend(legend_txt); adjustlegend(splt);    
    
    % Link axis limits of all the figures
    if isfield(option, 'linkfig')
        if ~strcmp(option.linkfig,'') | ~isempty(option.linkfig)
            linkaxes(ax,option.linkfig)
            % LinkFigures(ax,option.linkfig)
        end
    end
end