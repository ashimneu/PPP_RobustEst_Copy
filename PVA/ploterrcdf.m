function [fig] = ploterrcdf(output,option)
% INPUT:    solvername choices: td, lts, raps
%           option struct should include:
%           one of the properties: TDn_i, LTSn_i, RAPSn_i, default = 1

fig    = [];    
figtag = "errcdf";
fignum = getfignum(figtag);
eb_smaller_plot = 0; % enable smaller zoomed plot
p = output.p;   % parameters variable
trustprior = 0; % toggle such that all priors are always used / not

% toggle visibility
vis.ls   = 'on'; 
vis.lts  = 'on';
vis.raps = 'on';
vis.td   = 'on';

% colors
clr.ls = 'm';
clr.lts = 'r';
clr.raps = [0.5 0.5 0];
clr.td = 'c';

% curve markers
mrkr.ls = 'o';
mrkr.lts = '*';
mrkr.raps = 's';
mrkr.td = '^';

if p.eb_LTS  == 1 
    LTSn_i  = option.LTSn_i; % LTSn_i  - ith LTS estimation data  
    if LTSn_i == 1
        LTS_coverage = 0.75;
    elseif LTSn_i == 2
        LTS_coverage = 0.5;
    elseif LTSn_i  == 3
        LTS_coverage = 0.75;
    elseif LTSn_i  == 4
        LTS_coverage = 0.85;
    elseif LTSn_i  == 5
        LTS_coverage = 0.95;
    end
end 

%% ======================= ECEF Error CDF Curves ===========================
fig = figure(fignum); clf; hold on; grid on;
xlim([0.8e-1 0.5e1])
ylim([0 1.03])

legend_list = {};

% Least Squares Curve -------------------------------------------|
err_LS   = output.err_LS;
[f_ls,x_ls] = ecdf(err_LS);
h_ls = semilogx(x_ls,f_ls); hold on;
h_ls.MarkerSize   = 4;
legend_list = insert_legend('LS',legend_list,-1); 
set(h_ls,'LineWidth',1)
markeridx = 1:300:numel(h_ls.XData);
set(h_ls,'MarkerIndices',markeridx);

% tail stopper
stopper_line_width1 = 1.5;
stop_w1  = .02; % (1/2)width
stop_ls1 = line([x_ls(end) x_ls(end)],[1+stop_w1 1-stop_w1],...
    'HandleVisibility','off','LineWidth',stopper_line_width1);

% For Least squares line in main plot & small plot
set(h_ls,'Color',clr.ls); % line color
set(stop_ls1,'Color',clr.ls); % tail stopper color
set(h_ls,'Marker',mrkr.ls); % markers symbol
set(h_ls,'Visible',vis.ls); % cdf curve visibility
set(h_ls,'HandleVisibility',vis.ls); % legend visibility 
%-------------------------------------------------------------------|

if p.eb_LTS
    if trustprior, err_LTS = output.err_LTS2; else, err_LTS = output.err_LTS; end 
    [f_lts,x_lts] = ecdf(err_LTS(LTSn_i,:));
    h_lts = semilogx(x_lts,f_lts); hold on;
    h_lts.MarkerSize  = 5;
    legend_list = insert_legend(['LTS (g = floor(',num2str(LTS_coverage),'p))'],legend_list,-1);
    set(h_lts,'LineWidth',1)
    markeridx = 1:300:numel(h_ls.XData);
    set(h_lts,'MarkerIndices',markeridx);
    stop_lts1 = line([x_lts(end) x_lts(end)],[1+stop_w1 1-stop_w1],...
        'HandleVisibility','off','LineWidth',stopper_line_width1);
    set(h_lts,'Color',clr.lts);
    set(stop_lts1,'Color',clr.lts);
    set(h_lts,'Marker',mrkr.lts);    
    set(h_lts,'Visible',vis.lts);
    set(h_lts,'HandleVisibility',vis.lts); 
end

if p.eb_RAPS
    err_RAPS = output.err_RAPS;
    RAPSn_i  = option.RAPSn_i;
    [f_raps,x_raps] = ecdf(err_RAPS(RAPSn_i,:));
    h_raps = semilogx(x_raps,f_raps); hold on;
    h_raps.MarkerSize = 7;
    legend_list = insert_legend(['RAPS (',num2str(output.p.RAPSEps(RAPSn_i,1)),'m, ',num2str(output.p.RAPSEps(RAPSn_i,2)),'%)'],legend_list,-1);
    set(h_raps,'LineWidth',1)
    markeridx = 1:300:numel(h_ls.XData);
    set(h_raps,'MarkerIndices',markeridx);
    stop_raps1 = line([x_raps(end) x_raps(end)],[1+stop_w1 1-stop_w1],...
        'HandleVisibility','off','LineWidth',stopper_line_width1);
    set(h_raps,'Color',clr.raps);
    set(stop_raps1,'Color',clr.raps);
    set(h_raps,'Marker',mrkr.raps);    
    set(h_raps,'Visible',vis.raps)
    set(h_raps,'HandleVisibility',vis.raps)
end

if p.eb_TD
    err_TD = output.err_TD;
    TDn_i  = option.TDn_i;
    [f_td,x_td] = ecdf(err_TD(TDn_i,:));
    h_td = semilogx(x_td,f_td); hold on;
    h_td.MarkerSize   = 9;
    legend_list = insert_legend(['TD (\lambda = ',num2str(output.TDLambda(TDn_i)),')'],legend_list,-1);
    set(h_td,'LineWidth',1)
    markeridx = 1:300:numel(h_ls.XData);
    set(h_td,'MarkerIndices',markeridx);
    stop_td1   = line([x_td(end) x_td(end)],[1+stop_w1 1-stop_w1],...
        'HandleVisibility','off','LineWidth',stopper_line_width1);
    set(h_td,'Color',clr.td);
    set(stop_td1,'Color',clr.td);
    set(h_td,'Marker',mrkr.td);     
    set(h_td,'Visible',vis.td);
    set(h_td,'HandleVisibility',vis.td);
end

if p.eb_MShb
    if trustprior, err_MShb = output.err_MShb2; else, err_MShb = output.err_MShb; end 
    [f_mshb,x_mshb] = ecdf(err_MShb(MShbn_i,:));
    h_mshb = semilogx(x_mshb,f_mshb); hold on;
    h_mshb.MarkerSize  = 5;
    legend_list = insert_legend(['MER-Huber (c = ',num2str(output.MShbConst(MShbn_i)),')'],legend_list,-1);
    set(h_mshb,'LineWidth',1)
    markeridx = 1:300:numel(h_ls.XData);
    set(h_mshb,'MarkerIndices',markeridx);
    stop_mshb1 = line([x_mshb(end) x_mshb(end)],[1+stop_w1 1-stop_w1],...
        'HandleVisibility','off','LineWidth',stopper_line_width1);
    set(h_mshb,'Color',clr.mshb);
    set(stop_mshb1,'Color',clr.mshb);
    set(h_mshb,'Marker',mrkr.mshb);    
    set(h_mshb,'Visible',vis.mshb);
    set(h_mshb,'HandleVisibility',vis.mshb); 
end

if p.eb_MStk
    if trustprior, err_MStk = output.err_MStk2; else, err_MStk = output.err_MStk; end 
    [f_mstk,x_mstk] = ecdf(err_MStk(MStkn_i,:));
    h_mstk = semilogx(x_mstk,f_mstk); hold on;
    h_mstk.MarkerSize  = 5;
    legend_list = insert_legend(['MER-Tukey (c = ',num2str(output.MStkConst(MStkn_i)),')'],legend_list,-1);
    set(h_mstk,'LineWidth',1)
    markeridx = 1:300:numel(h_ls.XData);
    set(h_mstk,'MarkerIndices',markeridx);
    stop_mstk1 = line([x_mstk(end) x_mstk(end)],[1+stop_w1 1-stop_w1],...
        'HandleVisibility','off','LineWidth',stopper_line_width1);
    set(h_mstk,'Color',clr.mstk);
    set(stop_mstk1,'Color',clr.mstk);
    set(h_mstk,'Marker',mrkr.mstk);    
    set(h_mstk,'Visible',vis.mstk);
    set(h_mstk,'HandleVisibility',vis.mstk); 
end

if p.eb_LSS
    if trustprior, err_LSS = output.err_LSS2; else, err_LSS = output.err_LSS; end 
    [f_lss,x_lss] = ecdf(err_LSS(LSSn_i,:));
    h_lss = semilogx(x_lss,f_lss); hold on;
    h_lss.MarkerSize  = 5;
    legend_list = insert_legend(['LSS (\nu = ',num2str(output.LSSLambda(LSSn_i)),')'],legend_list,-1);
    set(h_lss,'LineWidth',1)
    markeridx = 1:300:numel(h_ls.XData);
    set(h_lss,'MarkerIndices',markeridx);
    stop_lss1 = line([x_lss(end) x_lss(end)],[1+stop_w1 1-stop_w1],...
        'HandleVisibility','off','LineWidth',stopper_line_width1);
    set(h_lss,'Color',clr.lss);
    set(stop_lss1,'Color',clr.lss);
    set(h_lss,'Marker',mrkr.lss);    
    set(h_lss,'Visible',vis.lss);
    set(h_lss,'HandleVisibility',vis.lss); 
end

Legend = legend(legend_list);
set(Legend,'Location','northwest'); set(Legend,'Interpreter','tex');

xlabel('Positioning Error (m)');
ylabel('Cumulative Probability');



if eb_smaller_plot
% ----------------------------Smaller Plot---------------------------------
%--------------------------------------------------------------------------

% rectangle box
rec_x = 1.38;       rec_width  = 25 - rec_x;
rec_y = 0.888;      rec_height = 1.01 - rec_y;
rec = rectangle('Position',[rec_x rec_y rec_width rec_height],'Curvature',0);
set(rec,'LineStyle',':');

% line connecting rectangle & smaller plot
line_left = line([1.38 1.90],[0.887 0.09],'HandleVisibility','off',...
    'LineWidth',2,'Color',[1 0.5 0],'LineStyle',':'); % left line [top_x bottom_x],[top_y bottom_y]
line_right = line([25 34],[1.01 0.84],'HandleVisibility','off',...
    'LineWidth',2,'Color',[1 0.5 0],'LineStyle',':'); % right line [top_x bottom_x],[top_y bottom_y]


% Axis
ax2 = axes('Position',[.53 .18 .37 .60]); hold on; box on; grid on;
axis([1.38 25 0.888 1.01]) % x-axis y-axis limits

% reference curve x values on which the smaller plot is based on
xref = x_ls; 
x1 = 1.3; 
x2 = 24.18;

% Least Squares Curve -------------------------------------------|
[~,idx1_ls] = min(abs(xref-x1)); 
[~,idx2_ls] = min(abs(xref-x2));
x_ls2 = x_ls(idx1_ls:1:idx2_ls); 
f_ls2 = f_ls(idx1_ls:1:idx2_ls);
h_ls2 = semilogx(x_ls2,f_ls2);

% adjust marker position
markeridx2 = 1:100:numel(h_ls2.XData); 
set(h_ls2,'MarkerIndices',markeridx2);

% line width
line_width2 = 0.7;
set(h_ls2,'LineWidth',line_width2)

% tail stopper
stopper_line_width22 = 1.5;
stop_w2    = .004; % (1/2)width
stop_ls2   = line([x_ls(end) x_ls(end)],[1+stop_w2 1-stop_w2],...
    'HandleVisibility','off','LineWidth',stopper_line_width22);

set(h_ls2,'Color',clr.ls); % line color
set(stop_ls2,'Color',clr.ls); % tail stopper color
set(h_ls2,'Marker',mrkr.ls); % markers symbol
set(h_ls2,'Visible',vis.ls); % cdf curve visibility
set(h_ls2,'HandleVisibility',vis.ls); % legend visibility 
%-------------------------------------------------------------------|

if p.eb_LTS
    [~,idx1_lts] = min(abs(xref-x1)); 
    [~,idx2_lts] = min(abs(xref-x2));
    x_lts2 = x_lts(idx1_lts:1:idx2_lts);  
    f_lts2 = f_lts(idx1_lts:1:idx2_lts);
    h_lts2 = semilogx(x_lts2,f_lts2);
    stop_lts2  = line([x_lts(end) x_lts(end)],[1+stop_w2 1-stop_w2],...
        'HandleVisibility','off','LineWidth',stopper_line_width22);
    markeridx2 = 1:100:numel(h_lts2.XData); 
    set(h_lts2,'MarkerIndices',markeridx2);
    set(h_lts2,'LineWidth',line_width2)
    set(h_lts2,'Color',clr.lts);
    set(stop_lts2,'Color',clr.lts);
    set(h_lts2,'Marker',mrkr.lts);    
    set(h_lts2,'Visible',vis.lts);
    set(h_lts2,'HandleVisibility',vis.lts);
end

if p.eb_RAPS
    [~,idx1_raps] = min(abs(xref-x1)); 
    [~,idx2_raps] = min(abs(xref-x2));
    x_raps2 = x_raps(idx1_raps:1:idx2_raps); 
    f_raps2 = f_raps(idx1_raps:1:idx2_raps);
    h_raps2 = semilogx(x_raps2,f_raps2);
    stop_raps2 = line([x_raps(end) x_raps(end)],[1+stop_w2 1-stop_w2],...
        'HandleVisibility','off','LineWidth',stopper_line_width22);
    markeridx2 = 1:100:numel(h_raps2.XData); 
    set(h_raps2,'MarkerIndices',markeridx2);
    set(h_raps2,'LineWidth',line_width2)
    set(h_raps2,'Color',clr.raps);
    set(stop_raps2,'Color',clr.raps);
    set(h_raps2,'Marker',mrkr.raps);    
    set(h_raps2,'Visible',vis.raps)
    set(h_raps2,'HandleVisibility',vis.raps)
end

if p.eb_TD
    [~,idx1_td] = min(abs(xref-x1)); 
    [~,idx2_td] = min(abs(xref-x2));
    x_td2 = x_td(idx1_td:1:idx2_td);       
    f_td2 = f_td(idx1_td:1:idx2_td);
    h_td2 = semilogx(x_td2,f_td2);
    stop_td2 = line([x_td(end) x_td(end)],[1+stop_w2 1-stop_w2],...
        'HandleVisibility','off','LineWidth',stopper_line_width22);
    markeridx2 = 1:100:numel(h_td2.XData);
    set(h_td2,'MarkerIndices',markeridx2);
    set(h_td2,'LineWidth',line_width2)
    set(h_td2,'Color',clr.td);
    set(stop_td2,'Color',clr.td);
    set(h_td2,'Marker',mrkr.td);     
    set(h_td2,'Visible',vis.td);
    set(h_td2,'HandleVisibility',vis.td);
end

end

% % adjust cdf figure window position
% set(fig,'Position',[20 50 640 480]);


end

function [adjusted_data,cutoff_data,all_flags,flag_sum] = cutoff(raw_data,threshold,n)
% trim values in raw_data such that entries larger than threshold
% are set to threshold value

% input: n - number of times each solver need

flag_idx = [];
cutoff_data = nan(size(raw_data));
adjusted_data = raw_data;
if ~isinf(threshold)    
    flag_idx = find(raw_data > threshold);
    adjusted_data(flag_idx) = threshold;   
end

cutoff_data(flag_idx) = threshold;
sum_flag = sum(~isnan(cutoff_data));

if nargout == 3
    all_flags = flag_idx;
elseif nargout == 4
    all_flags = flag_idx;
    flag_sum = sum_flag;
end

end
