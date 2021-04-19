function option = plotErrorStd(output,option)
    option = plotErrorStdSolver(output,"kf",option);
    if output.p.eb_LTS, option = plotErrorStdSolver(output,"lts",option); end
    if output.p.eb_TD,  option = plotErrorStdSolver(output,"td",option); end
    if output.p.eb_RAPS, option = plotErrorStdSolver(output,"raps",option); end
    if output.p.eb_MShb, option = plotErrorStdSolver(output,"mshb",option); end
    if output.p.eb_MStk, option = plotErrorStdSolver(output,"mstk",option); end
    
end

function option = plotErrorStdSolver(output,solvername,option)
    
    ax = [];    
    % Frame Conversion
    wgs84  = wgs84Ellipsoid(); % meters
    if option.movingrover == 1
        N = size(Grdpos,2); % GPS epochs
        Grdpos = output.p.Grdpos.pos;  % Ground Truth pos in ECEF
    else
        N = output.NN; %4475;
        Grdpos = repmat(output.p.Grdpos.pos(:,1),1,N);  % Ground Truth pos in ECEF
    end      
    grd0   = ecef2lla(Grdpos(:,1)'); % [degrees] true initial pos
    [trueN,trueE,trueD] = ecef2ned(Grdpos(1,:),Grdpos(2,:),Grdpos(3,:),grd0(1),grd0(2),grd0(3),wgs84);
    
    % Origin of NED frame in ECEF frame, i.e. True Position at time = 0sec
    NED_origin = Grdpos(:,1);
    p1_lla  = ecef2lla(Grdpos(:,1)'); % [degrees] true initial pos
    Rot = zeros(3,3); % Rotation matrix from ECEF to NED frame
    [Rot(1,1),Rot(2,1),Rot(3,1)] = ecef2ned(NED_origin(1)+1,NED_origin(2),NED_origin(3),p1_lla(1),p1_lla(2),p1_lla(3),wgs84);
    [Rot(2,1),Rot(2,2),Rot(3,2)] = ecef2ned(NED_origin(1),NED_origin(2)+1,NED_origin(3),p1_lla(1),p1_lla(2),p1_lla(3),wgs84);
    [Rot(1,3),Rot(2,3),Rot(3,3)] = ecef2ned(NED_origin(1),NED_origin(2),NED_origin(3)+1,p1_lla(1),p1_lla(2),p1_lla(3),wgs84);
    
    pltTstart1 = 1;          % <--------------- plot time start     
    pltTend1   = output.NN;  % <--------------- plot time end
    TW1 = pltTstart1:pltTend1; % time window 1
    postTime  = output.gpst(TW1); % begins at 0
      
    
    if lower(solvername) == "kf"
        solvername = "kf";
        post_cov  = output.post_statecovar.LS;
        err_3d = output.err_LS;
    
    elseif lower(solvername) == "lts"
        LTSn_i   = option.LTSn_i;
        post_cov  = output.post_statecovar.LTS{LTSn_i};
        err_3d = output.err_LTS(LTSn_i,:);
            
    elseif lower(solvername) == "td"
        TDn_i   = option.TDn_i;
        post_cov  = output.post_statecovar.TD{TDn_i};
        err_3d = output.err_TD(TDn_i,:);
            
    elseif lower(solvername) == "raps"
        RAPSn_i   = option.RAPSn_i;
        post_cov  = output.post_statecovar.RAPS{RAPSn_i};
        err_3d = output.err_RAPS(RAPSn_i,:);
        
    elseif lower(solvername) == "mshb"
        MShbn_i   = option.MShbn_i;
        post_cov  = output.post_statecovar.MShb{MShbn_i};
        err_3d = output.err_MShb(MShbn_i,:);
        
    elseif lower(solvername) == "mstk"
        MStkn_i   = option.MStkn_i;
        post_cov  = output.post_statecovar.MStk{MStkn_i};
        err_3d = output.err_MStk(MStkn_i,:);        
    end
    

    % Convert state covariances from ECEF to NED frame
    pos_post_cov  = pagemtimes(pagemtimes(Rot,post_cov(1:3,1:3,:)),Rot');    

    pspx_cov = squeeze(pos_post_cov(1,1,:));   
    pspy_cov = squeeze(pos_post_cov(2,2,:));  
    pspz_cov = squeeze(pos_post_cov(3,3,:)); 
    
    sampled_std   = nanstd(err_3d).*ones(1,N);
    estimated_std = sqrt(pspx_cov + pspy_cov + pspz_cov);
    
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
% Estimation & Predicted Covar
%--------------------------------------------------------------------------
    figtag = strcat("estpred_std",solvername);
    fignum = getfignum(figtag);    
    figure(fignum); clf; hold on; grid on
    xlbl_txt = strcat('Receiver time using GPS second');
    legend_txt = {'computed for all epoch','prediction at each epoch'};
    
    scatter(postTime,sampled_std(TW1),'Marker', mrkr,'MarkerEdgeColor',ps_mrk_clr,'SizeData',sz);
    scatter(postTime,estimated_std(TW1),'Marker', mrkr,'MarkerEdgeColor',psstd_mrk_clr,'SizeData',sz);
    ylabel('std, m')
    xlim(xlimits); ylim(ylimits)
    xlabel(xlbl_txt)    
    option.ax = [ax; findall(gcf, 'type', 'axes')];
    title_txt = strcat('Std. of estimation error (',upper(solvername),')');
    sgtitle(title_txt)
    legend(legend_txt); 
%     adjustlegend(gcf);
end

