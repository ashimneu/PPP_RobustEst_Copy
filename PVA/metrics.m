function [] = metrics(output,option)
    
    eb_GDOP = false;
    
    % Do ECEF error & GDOP cutoff
    err_cutoff  = 3;
    GDOP_cutoff = 3;

    err_LS = output.err_LS;     GDOPLS   = output.GDOPLS; 	nsvLS     = output.nsvLS;  
    dnsvLS = output.dnsvLS;     npriorLS = output.npriorLS; dnpriorLS = output.dnpriorLS;
    [ctf.err_LS,ctf.idx_err_LS,~,ctf.fsum_err_LS] = cutoff(err_LS,err_cutoff);
    [ctf.GDOPLS,ctf.idx_GDOPLS,~,ctf.fsum_GDOPLS] = cutoff(GDOPLS,GDOP_cutoff);
    
    if output.p.eb_TD
        TDn_i   = option.TDn_i;
        err_TD = output.err_TD; GDOPTD   = output.GDOPTD;   nsvTD     = output.nsvTD;  
        dnsvTD = output.dnsvTD; npriorTD = output.npriorTD; dnpriorTD = output.dnpriorTD;
        [ctf.err_TD,ctf.idx_err_TD,~,ctf.fsum_err_TD] = cutoff(err_TD(TDn_i,:),err_cutoff);
        [ctf.GDOPTD,ctf.idx_GDOPTD,~,ctf.fsum_GDOPTD] = cutoff(GDOPTD(TDn_i,:),GDOP_cutoff);
    end

    if output.p.eb_RAPS
        RAPSn_i = option.RAPSn_i;
        err_RAPS = output.err_RAPS; GDOPRAPS   = output.GDOPRAPS;   nsvRAPS     = output.nsvRAPS; 
        dnsvRAPS = output.dnsvRAPS; npriorRAPS = output.npriorRAPS; dnpriorRAPS = output.dnpriorRAPS;
        [ctf.err_RAPS,ctf.idx_err_RAPS,~,ctf.fsum_err_RAPS] = cutoff(err_RAPS(RAPSn_i,:),err_cutoff);
        [ctf.GDOPRAPS,ctf.idx_GDOPRAPS,~,ctf.fsum_GDOPRAPS] = cutoff(GDOPRAPS(RAPSn_i,:),GDOP_cutoff);
    end

    if output.p.eb_LTS
        LTSn_i  = option.LTSn_i;
        err_LTS  = output.err_LTS; GDOPLTS   = output.GDOPLTS;   nsvLTS     = output.nsvLTS;   
        dnsvLTS  = output.dnsvLTS; npriorLTS = output.npriorLTS; dnpriorLTS = output.dnpriorLTS;
        [ctf.err_LTS,ctf.idx_err_LTS,~,ctf.fsum_err_LTS] = cutoff(err_LTS(LTSn_i,:),err_cutoff);
        [ctf.GDOPLTS,ctf.idx_GDOPLTS,~,ctf.fsum_GDOPLTS] = cutoff(GDOPLTS(LTSn_i,:),GDOP_cutoff);
    end
    
    if output.p.eb_MShb
        MShbn_i = option.MShbn_i;
        err_MShb  = output.err_MShb; GDOPMShb   = output.GDOPMShb;   nsvMShb     = output.nsvMShb;   
        dnsvMShb  = output.dnsvMShb; npriorMShb = output.npriorMShb; dnpriorMShb = output.dnpriorMShb;
        [ctf.err_MShb,ctf.idx_err_MShb,~,ctf.fsum_err_MShb] = cutoff(err_MShb(MShbn_i,:),err_cutoff);
        [ctf.GDOPMShb,ctf.idx_GDOPMShb,~,ctf.fsum_GDOPMShb] = cutoff(GDOPMShb(MShbn_i,:),GDOP_cutoff);
    end
    
    if output.p.eb_MStk
        MStkn_i = option.MStkn_i;
        err_MStk  = output.err_MStk; GDOPMStk   = output.GDOPMStk;   nsvMStk     = output.nsvMStk;   
        dnsvMStk  = output.dnsvMStk; npriorMStk = output.npriorMStk; dnpriorMStk = output.dnpriorMStk;
        [ctf.err_MStk,ctf.idx_err_MStk,~,ctf.fsum_err_MStk] = cutoff(err_MStk(MStkn_i,:),err_cutoff);
        [ctf.GDOPMStk,ctf.idx_GDOPMStk,~,ctf.fsum_GDOPMStk] = cutoff(GDOPMStk(MStkn_i,:),GDOP_cutoff);
    end
    
    if output.p.eb_LSS
        LSSn_i  = option.LSSn_i;
        err_LSS  = output.err_LSS; GDOPLSS   = output.GDOPLSS;   nsvLSS     = output.nsvLSS;   
        dnsvLSS  = output.dnsvLSS; npriorLSS = output.npriorLSS; dnpriorLSS = output.dnpriorLSS;
        [ctf.err_LSS,ctf.idx_err_LSS,~,ctf.fsum_err_LSS] = cutoff(err_LSS(LSSn_i,:),err_cutoff);
        [ctf.GDOPLSS,ctf.idx_GDOPLSS,~,ctf.fsum_GDOPLSS] = cutoff(GDOPLSS(LSSn_i,:),GDOP_cutoff);
    end
    
    N = size(err_LS(1,:),2);
    fprintf('\nTotal epochs = %1.0f \n',N)

    % Error Standard Deviation
    std_LS = nanstd(err_LS(1,:));
    if output.p.eb_LTS,  std_LTS  = nanstd(err_LTS(LTSn_i,:));   end
    if output.p.eb_TD,   std_TD   = nanstd(err_TD(TDn_i,:));     end
    if output.p.eb_RAPS, std_RAPS = nanstd(err_RAPS(RAPSn_i,:)); end
    if output.p.eb_MShb, std_MShb = nanstd(err_MShb(MShbn_i,:)); end
    if output.p.eb_MStk, std_MStk = nanstd(err_MStk(MStkn_i,:)); end
    if output.p.eb_LSS,  std_LSS  = nanstd(err_LSS(LSSn_i,:));   end
    
    fprintf('\nError standard deviation (meters) \n')
    fprintf('LS   = %5.2f \n',std_LS)
    if output.p.eb_LTS,  fprintf('LTS  = %5.2f \n',std_LTS),    end
    if output.p.eb_TD,   fprintf('TD   = %5.2f \n',std_TD),     end
    if output.p.eb_RAPS, fprintf('RAPS = %5.2f \n',std_RAPS),   end
    if output.p.eb_MShb, fprintf('MShb = %5.2f \n',std_MShb),   end
    if output.p.eb_MStk, fprintf('MStk = %5.2f \n',std_MStk),   end
    if output.p.eb_LSS,  fprintf('LSS  = %5.2f \n',std_LSS),    end
    
    % Largest error
    fprintf('\nLargest Positioning error (meters) \n')
    fprintf('LS   = %5.2f \n',max(err_LS))
    if output.p.eb_LTS,  fprintf('LTS  = %5.2f \n',max(err_LTS(LTSn_i,:))),     end
    if output.p.eb_TD,   fprintf('TD   = %5.2f \n',max(err_TD(TDn_i,:))),       end
    if output.p.eb_RAPS, fprintf('RAPS = %5.2f \n',max(err_RAPS(RAPSn_i,:))),   end
    if output.p.eb_MShb, fprintf('MShb = %5.2f \n', max(err_MShb(MShbn_i,:))),end
    if output.p.eb_MStk, fprintf('MStk = %5.2f \n', max(err_MStk(MStkn_i,:))),end
    if output.p.eb_LSS,  fprintf('LSS  = %5.2f \n', max(err_LSS(LSSn_i,:))),  end
    
    % Error Prediction accuracy (probability of meeting predicted error)
    sfig = 2; % significant figure
    pred_acc_LS = round(sum(err_LS <= GDOPLS)/N,sfig)*100;
    if output.p.eb_LTS,  pred_acc_LTS  = round(sum(err_LTS(LTSn_i,:) <= GDOPLTS(LTSn_i,:))/N,sfig)*100;     end
    if output.p.eb_TD,   pred_acc_TD   = round(sum(err_TD(TDn_i,:) <= GDOPTD(TDn_i,:))/N,sfig)*100;         end
    if output.p.eb_RAPS, pred_acc_RAPS = round(sum(err_RAPS(RAPSn_i,:) <= GDOPRAPS(RAPSn_i,:))/N,sfig)*100; end
    if output.p.eb_MShb, pred_acc_MShb = round(sum(err_MShb(MShbn_i,:) <= GDOPMShb(MShbn_i,:))/N,sfig)*100; end
    if output.p.eb_MStk, pred_acc_MStk = round(sum(err_MStk(MStkn_i,:) <= GDOPMStk(MStkn_i,:))/N,sfig)*100; end
    if output.p.eb_LSS,  pred_acc_LSS  = round(sum(err_LSS(LSSn_i,:)  <= GDOPLSS(LSSn_i,:))/N,sfig)*100;   end
    
    fprintf('\nFraction of epochs where \npos error were cutoff at %2.1f meters. \n',err_cutoff)
    fprintf('LS   = %3.1f %%\n', ctf.fsum_err_LS/N*100)
    if output.p.eb_LTS,  fprintf('LTS  = %3.1f %%\n', ctf.fsum_err_LTS/N*100),  end
    if output.p.eb_TD,   fprintf('TD   = %3.1f %%\n', ctf.fsum_err_TD/N*100),   end
    if output.p.eb_RAPS, fprintf('RAPS = %3.1f %%\n', ctf.fsum_err_RAPS/N*100), end
    if output.p.eb_MShb, fprintf('MShb = %3.1f %%\n', ctf.fsum_err_MShb/N*100), end
    if output.p.eb_MStk, fprintf('MStk = %3.1f %%\n', ctf.fsum_err_MStk/N*100), end
    if output.p.eb_LSS,  fprintf('LSS  = %3.1f %%\n', ctf.fsum_err_LSS/N*100),  end
    
    if eb_GDOP
        fprintf('\nError Prediction Accuracy \n')
        fprintf('LS   = %3.1f %%\n', pred_acc_LS)
        if output.p.eb_LTS,  fprintf('LTS  = %3.1f %%\n', pred_acc_LTS),    end
        if output.p.eb_TD,   fprintf('TD   = %3.1f %%\n', pred_acc_TD),     end
        if output.p.eb_RAPS, fprintf('RAPS = %3.1f %%\n', pred_acc_RAPS),   end
        if output.p.eb_MShb, fprintf('MShb = %3.1f %%\n', pred_acc_MShb),   end
        if output.p.eb_MStk, fprintf('MStk = %3.1f %%\n', pred_acc_MStk),   end
        if output.p.eb_LSS,  fprintf('LSS  = %3.1f %%\n', pred_acc_LSS),    end

        fprintf('\nFraction of epochs for which  \nMSE_y were cutoff. \n')
        fprintf('LS   = %3.1f %%\n',ctf.fsum_GDOPLS/N*100)
        if output.p.eb_LTS,  fprintf('LTS  = %3.1f %%\n',ctf.fsum_GDOPLTS/N*100),   end
        if output.p.eb_TD,   fprintf('TD   = %3.1f %%\n',ctf.fsum_GDOPTD/N*100),    end
        if output.p.eb_RAPS, fprintf('RAPS = %3.1f %%\n',ctf.fsum_GDOPRAPS/N*100),  end
        if output.p.eb_MShb, fprintf('MShb = %3.1f %%\n',ctf.fsum_GDOPMShb/N*100),  end
        if output.p.eb_MStk, fprintf('MStk = %3.1f %%\n',ctf.fsum_GDOPMStk/N*100),  end
        if output.p.eb_LSS,  fprintf('LSS  = %3.1f %%\n',ctf.fsum_GDOPLSS/N*100),   end
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
