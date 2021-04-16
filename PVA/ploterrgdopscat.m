function option = ploterrgdopscat(output,solvername,option)
% INPUT:    solvername choices: td, lts, raps
%           option struct should include:
%           one of the properties: TDn_i, LTSn_i, RAPSn_i, default = 1

eb_GDOP = false; % enable GDOP scatter 
eb_nsv  = true;
eb_dnsv = true; % enable # discarded measurements


fig = []; ax = [];
p = output.p;
if solvername == "td",    if p.eb_TD ~= 1,   return, end, end
if solvername == "lts",   if p.eb_LTS ~= 1,  return, end, end
if solvername == "raps",  if p.eb_RAPS ~= 1, return, end, end
if solvername == "mshb",  if p.eb_MShb ~= 1, return, end, end
if solvername == "mstk",  if p.eb_MStk ~= 1, return, end, end
if solvername == "lss",  if p.eb_LSS ~= 1, return, end, end

solvername  = lower(solvername);
solvername2 = convertStringsToChars(upper(solvername));

%% indices of chosen parameters
if solvername == "td",    TDn_i   = option.TDn_i;   end     % TDn_i   - ith TD estimation data
if solvername == "lts",   LTSn_i  = option.LTSn_i;  end     % LTSn_i  - ith LTS estimation data
if solvername == "raps",  RAPSn_i = option.RAPSn_i; end     % RAPSn_i - ith RAPS estimation data
if solvername == "mshb",  MShbn_i = option.MShbn_i; end     % MShbn_i - ith MShb estimation data
if solvername == "mstk",  MStkn_i = option.MStkn_i; end     % MStkn_i - ith MStk estimation data
if solvername == "lss",  LSSn_i = option.LSSn_i; end     % LSSn_i - ith LSS estimation data
trustprior = 0;

if solvername == "lts"
    LTS_coverage = 0;
    if LTSn_i == 1
        LTS_coverage = 0.75;
    elseif LTSn_i == 2
        LTS_coverage = 0.5;
    end
end

%% Cutoff error & MSE_y vals beyond a threshold for better scatter plot
err_cutoff  = inf; % error cutoff threshold
GDOP_cutoff = 3; % MSE_y cutoff threshold

if solvername == "kf"
    err_LS = output.err_LS; GDOPLS = output.GDOPLS; nsvLS = output.nsvLS; dnsvLS = output.dnsvLS; npriorLS = output.npriorLS; dnpriorLS = output.dnpriorLS;
    [ctf.err_LS,ctf.idx_err_LS,~,ctf.fsum_err_LS] = cutoff(err_LS,err_cutoff); [ctf.GDOPLS,ctf.idx_GDOPLS,~,ctf.fsum_GDOPLS] = cutoff(GDOPLS,GDOP_cutoff);
    max_nsv = max(nsvLS,[],'all')+1;
end

if solvername == "lts"
    if trustprior, err_LTS = output.err_LTS2; GDOPLTS = output.GDOPLTS2; nsvLTS = output.nsvLTS2; dnsvLTS = output.dnsvLTS2; npriorLTS = output.npriorLTS2; dnpriorLTS = output.dnpriorLTS2;
    else err_LTS = output.err_LTS; GDOPLTS = output.GDOPLTS; nsvLTS = output.nsvLTS; dnsvLTS = output.dnsvLTS; npriorLTS = output.npriorLTS; dnpriorLTS = output.dnpriorLTS; end
    
    [ctf.err_LTS,ctf.idx_err_LTS,~,ctf.fsum_err_LTS] = cutoff(err_LTS(LTSn_i,:),err_cutoff); [ctf.GDOPLTS,ctf.idx_GDOPLTS,~,ctf.fsum_GDOPLTS] = cutoff(GDOPLTS(LTSn_i,:),GDOP_cutoff);
    max_nsv = max(nsvLTS,[],'all')+1;
end

if solvername == "raps"
    err_RAPS = output.err_RAPS; GDOPRAPS = output.GDOPRAPS; nsvRAPS = output.nsvRAPS; dnsvRAPS = output.dnsvRAPS;  npriorRAPS = output.npriorRAPS; dnpriorRAPS = output.dnpriorRAPS;
    [ctf.err_RAPS,ctf.idx_err_RAPS,~,ctf.fsum_err_RAPS] = cutoff(err_RAPS(RAPSn_i,:),err_cutoff); [ctf.GDOPRAPS,ctf.idx_GDOPRAPS,~,ctf.fsum_GDOPRAPS] = cutoff(GDOPRAPS(RAPSn_i,:),GDOP_cutoff);    
    max_nsv = max(nsvRAPS,[],'all')+1;
end

if solvername == "td"
    err_TD = output.err_TD; GDOPTD = output.GDOPTD; nsvTD = output.nsvTD; dnsvTD = output.dnsvTD; npriorTD = output.npriorTD; dnpriorTD = output.dnpriorTD;
    [ctf.err_TD,ctf.idx_err_TD,~,ctf.fsum_err_TD] = cutoff(err_TD(TDn_i,:),err_cutoff); [ctf.GDOPTD,ctf.idx_GDOPTD,~,ctf.fsum_GDOPTD] = cutoff(GDOPTD(TDn_i,:),GDOP_cutoff);
    max_nsv = max(nsvTD,[],'all')+1;
end

if solvername == "mshb"
    err_MShb = output.err_MShb; GDOPMShb = output.GDOPMShb; nsvMShb = output.nsvMShb; dnsvMShb = output.dnsvMShb;  npriorMShb = output.npriorMShb; dnpriorMShb = output.dnpriorMShb;
    [ctf.err_MShb,ctf.idx_err_MShb,~,ctf.fsum_err_MShb] = cutoff(err_MShb(MShbn_i,:),err_cutoff); [ctf.GDOPMShb,ctf.idx_GDOPMShb,~,ctf.fsum_GDOPMShb] = cutoff(GDOPMShb(MShbn_i,:),GDOP_cutoff);    
    max_nsv = max(nsvMShb,[],'all')+1;
end

if solvername == "mstk"
    err_MStk = output.err_MStk; GDOPMStk = output.GDOPMStk; nsvMStk = output.nsvMStk; dnsvMStk = output.dnsvMStk;  npriorMStk = output.npriorMStk; dnpriorMStk = output.dnpriorMStk;
    [ctf.err_MStk,ctf.idx_err_MStk,~,ctf.fsum_err_MStk] = cutoff(err_MStk(MStkn_i,:),err_cutoff); [ctf.GDOPMStk,ctf.idx_GDOPMStk,~,ctf.fsum_GDOPMStk] = cutoff(GDOPMStk(MStkn_i,:),GDOP_cutoff);    
    max_nsv = max(nsvMStk,[],'all')+1;
end

if solvername == "lss"
    err_LSS = output.err_LSS; GDOPLSS = output.GDOPLSS; nsvLSS = output.nsvLSS; dnsvLSS = output.dnsvLSS;  npriorLSS = output.npriorLSS; dnpriorLSS = output.dnpriorLSS;
    [ctf.err_LSS,ctf.idx_err_LSS,~,ctf.fsum_err_LSS] = cutoff(err_LSS(LSSn_i,:),err_cutoff); [ctf.GDOPLSS,ctf.idx_GDOPLSS,~,ctf.fsum_GDOPLSS] = cutoff(GDOPLSS(LSSn_i,:),GDOP_cutoff);    
    max_nsv = max(nsvLSS,[],'all')+1;
end

% markers
sz = 20; %10; % size of the markers
nsvMarker      = '.';           dnsvMarker     = 's';
npriorMarker   = '.';           dnpriorMarker  = '.';
GDOPMarker     = 'r.';          errMarker      = '.';
ctf_errMarker  = '^';           ctf_GDOPMarker = 'o';

% marker colors
nsv_clr      = [0 0.5 0.99];     % cyanish
dnsv_clr     = [0 0.5 0];        % green
nprior_clr   = [0.5 0 0];        % brown
dnprior_clr  = [0 0.5 0];        % green
err_mrkr_clr = [0.55 0.55 0.55]; % gray
idx_err_clr  = [0.55 0.55 0.55]; % gray
idx_GDOP_clr = [1 0 0];          % 

% Axis properties
Lyaxis_clr   = [0 0 0];          % black
Ryaxis_clr   = [0 0.5 0];        % green
xtickformat = 'HH:MM';
tickangle   = 0;
gpst        = output.gpst;       % x axis values
xlimits     = [gpst(1) gpst(end)];
ylimit_l    = [0 err_cutoff+0.1];
ylimit_r    = [0 inf];
Lylbl_txt   = 'MSE_y & Positioning error';
Rylbl_txt 	= 'No. of measurements';
xlbl_txt    = 'Receiver time using GPS second'; %'local time (hr)';

%% Figure
figtag = strcat("errgdopnsv_",solvername);
fignum = getfignum(figtag);
fig    = figure(fignum); clf; hold on; grid on
legend_list = {};
yyaxis left; ylim(ylimit_l); ylabel(Lylbl_txt)

switch solvername
case "kf"    
    title({'KF Positioning error',['Measurement Update using ', solvername2 ]},'Interpreter','Latex')    
    if eb_GDOP
        sc6 = scatter(gpst,ctf.idx_GDOPLS,sz,ctf_GDOPMarker,'HandleVisibility','off');
        sc2 = scatter(gpst,ctf.GDOPLS,sz,GDOPMarker);         %#ok<*NASGU>
    end
    sc5 = scatter(gpst,ctf.idx_err_LS,sz,ctf_errMarker,'HandleVisibility','off');
    sc3 = scatter(gpst,ctf.err_LS,sz,errMarker);
    yyaxis right; ylim(ylimit_r)
    if eb_nsv,  sc1 = scatter(gpst,nsvLS,sz,nsvMarker);     end
    if eb_dnsv, sc4 = scatter(gpst,dnsvLS,sz,dnsvMarker);   end
    
case "lts"
    title({'ECEF Positioning error',['Measurement Update using ', solvername2,...
         ' (g = floor ',num2str(LTS_coverage),'p)']},'Interpreter','Latex')
    if eb_GDOP 
        sc6 = scatter(gpst,ctf.idx_GDOPLTS,sz,ctf_GDOPMarker,'HandleVisibility','off');
        sc2 = scatter(gpst,ctf.GDOPLTS,sz,GDOPMarker);        
    end
    sc5 = scatter(gpst,ctf.idx_err_LTS,sz,ctf_errMarker,'HandleVisibility','off');
    sc3 = scatter(gpst,ctf.err_LTS,sz,errMarker);
    yyaxis right; ylim(ylimit_r)
    if eb_nsv,  sc1 = scatter(gpst,nsvLTS(LTSn_i,:),sz,nsvMarker);      end
    if eb_dnsv, sc4 = scatter(gpst,dnsvLTS(LTSn_i,:),sz,dnsvMarker);    end
    
case "td"
    title({'ECEF Positioning error',['Measurement Update using ', solvername2,...
         ' ($\lambda$ = ',num2str(output.TDLambda(TDn_i)),')']},'Interpreter','Latex')    
    if eb_GDOP
        sc6 = scatter(gpst,ctf.idx_GDOPTD,sz,ctf_GDOPMarker,'HandleVisibility','off');
        sc2 = scatter(gpst,ctf.GDOPTD,sz,GDOPMarker);
    end
    sc5 = scatter(gpst,ctf.idx_err_TD,sz,ctf_errMarker,'HandleVisibility','off');
    sc3 = scatter(gpst,ctf.err_TD,sz,errMarker);
    yyaxis right; ylim(ylimit_r)
    if eb_nsv,  sc1 = scatter(gpst,nsvTD(TDn_i,:),sz,nsvMarker);       end
    if eb_dnsv, sc4 = scatter(gpst,dnsvTD(TDn_i,:),sz,dnsvMarker);     end
    
case "raps" 
    title({'ECEF Positioning error',['Measurement Update using ', solvername2,...
         ' (',num2str(output.p.RAPSEps(RAPSn_i,1)),'m, ',num2str(output.p.RAPSEps(RAPSn_i,2)),'\%)']},'Interpreter','Latex')
    if eb_GDOP
        sc6 = scatter(gpst,ctf.idx_GDOPRAPS,sz,ctf_GDOPMarker,'HandleVisibility','off');
        sc2 = scatter(gpst,ctf.GDOPRAPS,sz,GDOPMarker);
    end
    sc5 = scatter(gpst,ctf.idx_err_RAPS,sz,ctf_errMarker,'HandleVisibility','off');
    sc3 = scatter(gpst,ctf.err_RAPS,sz,errMarker);
    yyaxis right; ylim(ylimit_r)
    if eb_nsv,  sc1 = scatter(gpst,nsvRAPS(RAPSn_i,:),sz,nsvMarker);   end
    if eb_dnsv, sc4 = scatter(gpst,dnsvRAPS(RAPSn_i,:),sz,dnsvMarker); end

    case "mshb" 
    title({'ECEF Positioning error',['Measurement Update using ', solvername2,...
         ' (',num2str(output.p.MShbConst(MShbn_i,1)),')']},'Interpreter','Latex')
    if eb_GDOP
        sc6 = scatter(gpst,ctf.idx_GDOPMShb,sz,ctf_GDOPMarker,'HandleVisibility','off');
        sc2 = scatter(gpst,ctf.GDOPMShb,sz,GDOPMarker);
    end
    sc5 = scatter(gpst,ctf.idx_err_MShb,sz,ctf_errMarker,'HandleVisibility','off');
    sc3 = scatter(gpst,ctf.err_MShb,sz,errMarker);
    yyaxis right; ylim(ylimit_r)
    if eb_nsv,  sc1 = scatter(gpst,nsvMShb(MShbn_i,:),sz,nsvMarker);   end
    if eb_dnsv, sc4 = scatter(gpst,dnsvMShb(MShbn_i,:),sz,dnsvMarker); end

case "mstk" 
    title({'ECEF Positioning error',['Measurement Update using ', solvername2,...
         ' (',num2str(output.p.MStkConst(MStkn_i,1)),')']},'Interpreter','Latex')
    if eb_GDOP
        sc6 = scatter(gpst,ctf.idx_GDOPMStk,sz,ctf_GDOPMarker,'HandleVisibility','off');
        sc2 = scatter(gpst,ctf.GDOPMStk,sz,GDOPMarker);
    end
    sc5 = scatter(gpst,ctf.idx_err_MStk,sz,ctf_errMarker,'HandleVisibility','off');
    sc3 = scatter(gpst,ctf.err_MStk,sz,errMarker);
    yyaxis right; ylim(ylimit_r)
    if eb_nsv,  sc1 = scatter(gpst,nsvMStk(MStkn_i,:),sz,nsvMarker);   end
    if eb_dnsv, sc4 = scatter(gpst,dnsvMStk(MStkn_i,:),sz,dnsvMarker); end

case "lss" 
    title({'KF ECEF Positioning error',['Measurement Update using ', solvername2,...
         ' (',num2str(output.p.LSSLambda(LSSn_i,1)),')']},'Interpreter','Latex')
    if eb_GDOP
        sc6 = scatter(gpst,ctf.idx_GDOPLSS,sz,ctf_GDOPMarker,'HandleVisibility','off');
        sc2 = scatter(gpst,ctf.GDOPLSS,sz,GDOPMarker);
    end
    sc5 = scatter(gpst,ctf.idx_err_LSS,sz,ctf_errMarker,'HandleVisibility','off');
    sc3 = scatter(gpst,ctf.err_LSS,sz,errMarker);
    yyaxis right; ylim(ylimit_r)
    if eb_nsv,  sc1 = scatter(gpst,nsvLSS(LSSn_i,:),sz,nsvMarker);   end
    if eb_dnsv, sc4 = scatter(gpst,dnsvLSS(LSSn_i,:),sz,dnsvMarker); end
end

if eb_GDOP
    set(sc6,'CData',idx_GDOP_clr);
    legend_list = insert_legend('MSE_{y} (meter)',legend_list,-1);    
end

set(sc3,'CData',err_mrkr_clr);
set(sc5,'CData',idx_err_clr);
legend_list = insert_legend('Positioning error norm (meter)',legend_list,-1); 

if eb_nsv
    set(sc1,'CData',nsv_clr);
    legend_list = insert_legend('No. of used measurements',legend_list,-1);    
end
if eb_dnsv
    set(sc4,'CData',dnsv_clr);
    legend_list = insert_legend('No. of unused measurements',legend_list,-1);
end
    
ax = gca; 
set(ax.YAxis(1), 'Color',Lyaxis_clr);
set(ax.YAxis(2), 'Color',Ryaxis_clr);

ylabel(Rylbl_txt)
xlim(xlimits)
xlabel(xlbl_txt)

Legend = legend(legend_list);
set(Legend,'Orientation','Horizontal');
legend_position = [0.29,0.015,0.427,0.033];
Legend.Position = legend_position;

if ~isempty(ax) 
    if isfield(option,'ax')
        option.ax = [option.ax; ax];
    else
        option.ax = ax;
    end
end

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
