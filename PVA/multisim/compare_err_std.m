function [] = compare_err_std(outputcell,option,compar_param,compar_metric)

p = outputcell{1}.p;

if p.eb_outlier == 1
std_LS = [];
std_LTS = [];
std_RAPS = [];
std_TD = [];
std_MShb = [];
std_MStk = [];

maxerr_LS = [];
maxerr_LTS = [];
maxerr_RAPS = [];
maxerr_TD = [];
maxerr_MShb = [];
maxerr_MStk = [];

legend_list = {};

metrics_out_list = {}; 
simvarlist = [];
switch lower(compar_param)
    case "mean"
        simvarlist = p.outlierparam.meanlist;
        metrics_out_list = cell(numel(simvarlist),1);
        for sim_idx = 1:numel(simvarlist)
            output = outputcell{sim_idx};
            metrics_out_list{sim_idx} = metrics2(output,option);            
        end
        xlbl_txt = strcat('artificial added outlier mean (meters)');
    case "count"
        simvarlist = p.outlierparam.countlist;
        metrics_out_list = cell(numel(simvarlist),1);
        for sim_idx = 1:numel(simvarlist)
            output = outputcell{sim_idx};
            metrics_out_list{sim_idx} = metrics2(output,option);             
        end
        xlbl_txt = strcat('artifically added outlier count');
    case "width"
end

for j = 1:1:numel(metrics_out_list)
    metrics_out = metrics_out_list{j};
    std_LS = [std_LS metrics_out.std_LS];
    maxerr_LS = [maxerr_LS metrics_out.maxerr_LS];
    if p.eb_LTS  
        std_LTS = [std_LTS metrics_out.std_LTS];
        maxerr_LTS = [maxerr_LTS metrics_out.maxerr_LTS];
    end
    if p.eb_TD  
        std_TD = [std_TD metrics_out.std_TD];
        maxerr_TD = [maxerr_TD metrics_out.maxerr_TD];
    end
    if p.eb_RAPS
        std_RAPS = [std_RAPS metrics_out.std_RAPS]; 
        maxerr_RAPS = [maxerr_RAPS metrics_out.maxerr_RAPS];
    end
    if p.eb_MShb
        std_MShb = [std_MShb metrics_out.std_MShb]; 
        maxerr_MShb = [maxerr_MShb metrics_out.maxerr_MShb];
    end
    if p.eb_MStk 
        std_MStk = [std_MStk metrics_out.std_MStk];
        maxerr_MStk = [maxerr_MStk metrics_out.maxerr_MStk];
    end
end  

% switch lower(solvername)
%     case "kf"
%         err_std_solver = std_LS;
%         max_err_solver = maxerr_LS;
%     case "lts"
%         err_std_solver = std_LTS;
%         max_err_solver = maxerr_LTS;
%     case "raps"
%         err_std_solver = std_RAPS;
%         max_err_solver = maxerr_RAPS;
%     case "td"
%         err_std_solver = std_TD;
%         max_err_solver = maxerr_TD;
%     case "mshb"
%         err_std_solver = std_MShb;
%         max_err_solver = maxerr_MShb;
%     case "mstk"
%         err_std_solver = std_MStk;
%         max_err_solver = maxerr_MStk;
% end

figtag = strcat("compare_",compar_metric,"_",compar_param);
fignum = getfignum(figtag);    
figure(fignum); clf; hold on; grid on

%     if lower(compar_metric) == "std"
%         plot(simvarlist,err_std_solver,'.')
%     elseif lower(compar_metric) == "maxerr"
%         plot(simvarlist,max_err_solver,'.')
%     end
%    

% colors
clr.ls = 'm';
clr.lts = 'r';
clr.raps = [0.5 0.5 0];
clr.td = 'c';
clr.mshb = [0.9290 0.6940 0.1250];
clr.mstk = [0.3010 0.7450 0.9330];

% curve markers
mrkr.ls = 'o';
mrkr.lts = '*';
mrkr.raps = 's';
mrkr.td = '^';
mrkr.mshb = 'h';
mrkr.mstk = 'p';

if lower(compar_metric) == "std"
    ylbl_txt = strcat('positoning error std. (meters)');
%     plot(simvarlist,std_LS,'-')
%     if output.p.eb_LTS,  plot(simvarlist,std_LTS,'-');  end
%     if output.p.eb_TD,   plot(simvarlist,std_TD,'-');   end
%     if output.p.eb_RAPS, plot(simvarlist,std_RAPS,'-'); end
%     if output.p.eb_MShb, plot(simvarlist,std_MShb,'-'); end
%     if output.p.eb_MStk, plot(simvarlist,std_MStk,'-'); end
    plot(simvarlist,std_LS,'Marker',mrkr.ls)
    if output.p.eb_LTS,  plot(simvarlist,std_LTS,'Marker',mrkr.lts);  end
    if output.p.eb_TD,   plot(simvarlist,std_TD,'Marker',mrkr.td);   end
    if output.p.eb_RAPS, plot(simvarlist,std_RAPS,'Marker',mrkr.raps); end
    if output.p.eb_MShb, plot(simvarlist,std_MShb,'Marker',mrkr.mshb); end
    if output.p.eb_MStk, plot(simvarlist,std_MStk,'Marker',mrkr.mstk); end
    
elseif lower(compar_metric) == "maxerr"
    ylbl_txt = strcat('max positioning error (meters)');
    plot(simvarlist,maxerr_LS,'-')
    if output.p.eb_LTS,  plot(simvarlist,maxerr_LTS,'-');  end
    if output.p.eb_TD,   plot(simvarlist,maxerr_TD,'-');   end
    if output.p.eb_RAPS, plot(simvarlist,maxerr_RAPS,'-'); end
    if output.p.eb_MShb, plot(simvarlist,maxerr_MShb,'-'); end
    if output.p.eb_MStk, plot(simvarlist,maxerr_MStk,'-'); end
end

legend_list = insert_legend('LS',legend_list,-1); 
if output.p.eb_LTS,  legend_list = insert_legend('LTS',legend_list,-1);  end
if output.p.eb_TD,   legend_list = insert_legend('TD',legend_list,-1);   end
if output.p.eb_RAPS, legend_list = insert_legend('RAPS',legend_list,-1); end
if output.p.eb_MShb, legend_list = insert_legend('MShb',legend_list,-1); end
if output.p.eb_MStk, legend_list = insert_legend('MStk',legend_list,-1); end

Legend = legend(legend_list);
xlabel(xlbl_txt)
ylabel(ylbl_txt)

end
end

