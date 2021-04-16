function [] = metrics(output,option)
    
    eb_GDOP = false;
    
    % Do ECEF error & GDOP cutoff
    err_cutoff  = 3;
    GDOP_cutoff = 3;
    
    TDn_i   = option.TDn_i;
    LTSn_i  = option.LTSn_i;
    RAPSn_i = option.LTSn_i;
    trustprior = 0;

    err_LS = output.err_LS;     GDOPLS   = output.GDOPLS; 	nsvLS     = output.nsvLS;  
    dnsvLS = output.dnsvLS;     npriorLS = output.npriorLS; dnpriorLS = output.dnpriorLS;
    [ctf.err_LS,ctf.idx_err_LS,~,ctf.fsum_err_LS] = cutoff(err_LS,err_cutoff);
    [ctf.GDOPLS,ctf.idx_GDOPLS,~,ctf.fsum_GDOPLS] = cutoff(GDOPLS,GDOP_cutoff);
    
    if output.p.eb_TD
        err_TD = output.err_TD; GDOPTD   = output.GDOPTD;   nsvTD     = output.nsvTD;  
        dnsvTD = output.dnsvTD; npriorTD = output.npriorTD; dnpriorTD = output.dnpriorTD;
        [ctf.err_TD,ctf.idx_err_TD,~,ctf.fsum_err_TD] = cutoff(err_TD(TDn_i,:),err_cutoff);
        [ctf.GDOPTD,ctf.idx_GDOPTD,~,ctf.fsum_GDOPTD] = cutoff(GDOPTD(TDn_i,:),GDOP_cutoff);
    end

    if output.p.eb_RAPS
        err_RAPS = output.err_RAPS; GDOPRAPS   = output.GDOPRAPS;   nsvRAPS     = output.nsvRAPS; 
        dnsvRAPS = output.dnsvRAPS; npriorRAPS = output.npriorRAPS; dnpriorRAPS = output.dnpriorRAPS;
        [ctf.err_RAPS,ctf.idx_err_RAPS,~,ctf.fsum_err_RAPS] = cutoff(err_RAPS(RAPSn_i,:),err_cutoff);
        [ctf.GDOPRAPS,ctf.idx_GDOPRAPS,~,ctf.fsum_GDOPRAPS] = cutoff(GDOPRAPS(RAPSn_i,:),GDOP_cutoff);
    end

    if output.p.eb_LTS
    if trustprior == 1
        err_LTS  = output.err_LTS2; GDOPLTS   = output.GDOPLTS2;   nsvLTS     = output.nsvLTS2;
        dnsvLTS  = output.dnsvLTS2; npriorLTS = output.npriorLTS2; dnpriorLTS = output.dnpriorLTS2;
    else
        err_LTS  = output.err_LTS; GDOPLTS   = output.GDOPLTS;   nsvLTS     = output.nsvLTS;   
        dnsvLTS  = output.dnsvLTS; npriorLTS = output.npriorLTS; dnpriorLTS = output.dnpriorLTS;
    end
        [ctf.err_LTS,ctf.idx_err_LTS,~,ctf.fsum_err_LTS] = cutoff(err_LTS(LTSn_i,:),err_cutoff);
        [ctf.GDOPLTS,ctf.idx_GDOPLTS,~,ctf.fsum_GDOPLTS] = cutoff(GDOPLTS(LTSn_i,:),GDOP_cutoff);
    end

    
    N = size(err_LS(1,:),2);
    fprintf('\nTotal epochs = %1.0f \n',N)

    % Error Standard Deviation
                std_LS   = nanstd(err_LS(1,:));
    if output.p.eb_LTS,  std_LTS  = nanstd(err_LTS(1,:));        end
    if output.p.eb_TD,   std_TD   = nanstd(err_TD(TDn_i,:));     end
    if output.p.eb_RAPS, std_RAPS = nanstd(err_RAPS(RAPSn_i,:)); end

    fprintf('\nError standard deviation (meters) \n')
                fprintf('LS   = %5.2f \n',std_LS)
    if output.p.eb_LTS,  fprintf('LTS  = %5.2f \n',std_LTS),    end
    if output.p.eb_TD,   fprintf('TD   = %5.2f \n',std_TD),     end
    if output.p.eb_RAPS, fprintf('RAPS = %5.2f \n',std_RAPS),   end

    % Largest error
    fprintf('\nLargest Positioning error (meters) \n')
                fprintf('LS   = %5.2f \n',max(err_LS))
    if output.p.eb_LTS,  fprintf('LTS  = %5.2f \n',max(err_LTS(LTSn_i,:))),     end
    if output.p.eb_TD,   fprintf('TD   = %5.2f \n',max(err_TD(TDn_i,:))),       end
    if output.p.eb_RAPS, fprintf('RAPS = %5.2f \n',max(err_RAPS(RAPSn_i,:))),   end

    % Error Prediction accuracy (probability of meeting predicted error)
    sfig = 2; % significant figure
                pred_acc_LS = round(sum(err_LS <= GDOPLS)/N,sfig)*100;
    if output.p.eb_LTS,  pred_acc_LTS  = round(sum(err_LTS(LTSn_i,:) <= GDOPLTS(LTSn_i,:))/N,sfig)*100;     end
    if output.p.eb_TD,   pred_acc_TD   = round(sum(err_TD(TDn_i,:) <= GDOPTD(TDn_i,:))/N,sfig)*100;         end
    if output.p.eb_RAPS, pred_acc_RAPS = round(sum(err_RAPS(RAPSn_i,:) <= GDOPRAPS(RAPSn_i,:))/N,sfig)*100; end
    
    fprintf('\nFraction of epochs where \npos error were cutoff at %2.1f m. \n',err_cutoff)
                         fprintf('LS   = %3.1f %%\n', ctf.fsum_err_LS/N*100)
    if output.p.eb_LTS,  fprintf('LTS  = %3.1f %%\n', ctf.fsum_err_LTS/N*100),  end
    if output.p.eb_TD,   fprintf('TD   = %3.1f %%\n', ctf.fsum_err_TD/N*100),   end
    if output.p.eb_RAPS, fprintf('RAPS = %3.1f %%\n', ctf.fsum_err_RAPS/N*100), end
    
    if eb_GDOP
        fprintf('\nError Prediction Accuracy \n')
        fprintf('LS   = %3.1f %%\n', pred_acc_LS)
        if output.p.eb_LTS,  fprintf('LTS  = %3.1f %%\n', pred_acc_LTS),    end
        if output.p.eb_TD,   fprintf('TD   = %3.1f %%\n', pred_acc_TD),     end
        if output.p.eb_RAPS, fprintf('RAPS = %3.1f %%\n', pred_acc_RAPS),   end

        fprintf('\nFraction of epochs for which  \nMSE_y were cutoff. \n')
        fprintf('LS   = %3.1f %%\n',ctf.fsum_GDOPLS/N*100)
        if output.p.eb_LTS,  fprintf('LTS  = %3.1f %%\n',ctf.fsum_GDOPLTS/N*100),   end
        if output.p.eb_TD,   fprintf('TD   = %3.1f %%\n',ctf.fsum_GDOPTD/N*100),    end
        if output.p.eb_RAPS, fprintf('RAPS = %3.1f %%\n',ctf.fsum_GDOPRAPS/N*100),  end
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
