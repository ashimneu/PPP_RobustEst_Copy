if output.eb_LTS

    % Positioning Error
    fprintf('----------------- Positioning norm Error -------- \n')
    fprintf('standard deviation \n')
    LSstd = nanstd(output.err_LS(Timeline));
    LTSstd1 = nanstd(output.err_LTS(1,Timeline));
    LTSstd2 = nanstd(output.err_LTS(2,Timeline));
    LTSstd3 = nanstd(output.err_LTS(3,Timeline));
    LTSstd4 = nanstd(output.err_LTS(4,Timeline));
    fprintf('\nstd LS : %2.4f \n', LSstd)
    fprintf('std LTS regular selec. : %2.4f \n', LTSstd1)
    fprintf('std LTS, b = [b_(pse) & b_(dop)] : %2.4f \n', LTSstd2)
    fprintf('std LTS, b = [b_(pse) & b_(pse)] : %2.4f \n', LTSstd3)
    fprintf('std LTS, b = [b_(dop) & b_(dop)] : %2.4f \n', LTSstd4)

    fprintf('\nLargest Positioning error (meters) \n')
    LTS_maxerr_1 = max(output.err_LTS(1,Timeline));
    LTS_maxerr_2 = max(output.err_LTS(2,Timeline));
    LTS_maxerr_3 = max(output.err_LTS(3,Timeline));
    LTS_maxerr_4 = max(output.err_LTS(4,Timeline));
    fprintf('max err LS  = %5.2f \n',max(output.err_LS(startTime:endTime)))
    fprintf('max err LTS regular selec. : %2.4f \n', LTS_maxerr_1)
    fprintf('max err LTS, b = [b_(pse) & b_(dop)] : %2.4f \n', LTS_maxerr_2)
    fprintf('max err LTS, b = [b_(pse) & b_(pse)] : %2.4f \n', LTS_maxerr_3)
    fprintf('max err LTS, b = [b_(dop) & b_(dop)] : %2.4f \n', LTS_maxerr_4)

    err_cutoff = 3;
    [~,~,~,ctf_pcent_LS] = cutoff(output.err_LS(Timeline),err_cutoff);
    [~,~,~,ctf_pcent_LTS1] = cutoff(output.err_LTS(1,Timeline),err_cutoff);
    [~,~,~,ctf_pcent_LTS2] = cutoff(output.err_LTS(2,Timeline),err_cutoff);
    [~,~,~,ctf_pcent_LTS3] = cutoff(output.err_LTS(3,Timeline),err_cutoff);
    [~,~,~,ctf_pcent_LTS4] = cutoff(output.err_LTS(4,Timeline),err_cutoff);

    fprintf('\nFraction of epochs where \npos error were cutoff at %2.1f m. \n',err_cutoff)
    fprintf('LS   = %3.1f %%\n', ctf_pcent_LS/N*100)
    fprintf('LTS (regular) : %2.4f %%\n', ctf_pcent_LTS1/N*100)
    fprintf('LTS, b = [b_(pse) & b_(dop)] : %2.4f %%\n', ctf_pcent_LTS2/N*100)
    fprintf('LTS, b = [b_(pse) & b_(pse)] : %2.4f %%\n', ctf_pcent_LTS3/N*100)
    fprintf('LTS, b = [b_(dop) & b_(dop)] : %2.4f %%\n', ctf_pcent_LTS4/N*100)

    % Velocity Error
    fprintf('----------------- Velocity Norm Error -------- \n')
    verr_LS = output.post_state.LS{1};
    fprintf('standard deviation \n')
    LSstd = nanstd(output.err_LS(Timeline));
    LTSstd1 = nanstd(output.err_LTS(1,Timeline));
    LTSstd2 = nanstd(output.err_LTS(2,Timeline));
    LTSstd3 = nanstd(output.err_LTS(3,Timeline));
    LTSstd4 = nanstd(output.err_LTS(4,Timeline));
    fprintf('\nstd LS : %2.4f \n', LSstd)
    fprintf('std LTS regular selec. : %2.4f \n', LTSstd1)
    fprintf('std LTS, b = [b_(pse) & b_(dop)] : %2.4f \n', LTSstd2)
    fprintf('std LTS, b = [b_(pse) & b_(pse)] : %2.4f \n', LTSstd3)
    fprintf('std LTS, b = [b_(dop) & b_(dop)] : %2.4f \n', LTSstd4)

    fprintf('\nLargest Positioning error (meters) \n')
    LTS_maxerr_1 = max(output.err_LTS(1,Timeline));
    LTS_maxerr_2 = max(output.err_LTS(2,Timeline));
    LTS_maxerr_3 = max(output.err_LTS(3,Timeline));
    LTS_maxerr_4 = max(output.err_LTS(4,Timeline));
    fprintf('max err LS  = %5.2f \n',max(output.err_LS(startTime:endTime)))
    fprintf('max err LTS regular selec. : %2.4f \n', LTS_maxerr_1)
    fprintf('max err LTS, b = [b_(pse) & b_(dop)] : %2.4f \n', LTS_maxerr_2)
    fprintf('max err LTS, b = [b_(pse) & b_(pse)] : %2.4f \n', LTS_maxerr_3)
    fprintf('max err LTS, b = [b_(dop) & b_(dop)] : %2.4f \n', LTS_maxerr_4)

    err_cutoff = 3;
    [~,~,~,ctf_pcent_LS] = cutoff(output.err_LS(Timeline),err_cutoff);
    [~,~,~,ctf_pcent_LTS1] = cutoff(output.err_LTS(1,Timeline),err_cutoff);
    [~,~,~,ctf_pcent_LTS2] = cutoff(output.err_LTS(2,Timeline),err_cutoff);
    [~,~,~,ctf_pcent_LTS3] = cutoff(output.err_LTS(3,Timeline),err_cutoff);
    [~,~,~,ctf_pcent_LTS4] = cutoff(output.err_LTS(4,Timeline),err_cutoff);

    fprintf('\nFraction of epochs where \npos error were cutoff at %2.1f m. \n',err_cutoff)
    fprintf('LS   = %3.1f %%\n', ctf_pcent_LS/N*100)
    fprintf('LTS (regular) : %2.4f %%\n', ctf_pcent_LTS1/N*100)
    fprintf('LTS, b = [b_(pse) & b_(dop)] : %2.4f %%\n', ctf_pcent_LTS2/N*100)
    fprintf('LTS, b = [b_(pse) & b_(pse)] : %2.4f %%\n', ctf_pcent_LTS3/N*100)
    fprintf('LTS, b = [b_(dop) & b_(dop)] : %2.4f %%\n', ctf_pcent_LTS4/N*100)


    % nanmax(output.err_LTS(opt.LTSn_i,startTime:endTime))
    % nanstd(output.err(1:4300))
    % nanstd(output.vnorm(1:4300))

    err_diff = output.err_LTS(opt.LTSn_i,:) - output.err;
    pos_err_diff = err_diff >= 0.1;
    % count = sum(pos_err_diff)
    [max_res1,idx_max1] = maxk(err_diff,10);
    % max_res1 
    % idx_max1

    % figure(20); clf; hold on; grid on
    % plot(idx_max1,max_res1,'r*')
    % xlim([0 4475])

    byLTS = output.byLTS{opt.LTSn_i};
    byLTS = byLTS(any(byLTS,2),:);
end