function outputcell = compute_gnss_ecef_multisim(p,eph,obs)

outputcell = {};
if p.multisim == 0
    output = compute_gnss_ecef(p,eph,obs);
    outputcell{1} = output;
else
    switch lower(p.multisim_outliervar) 
        case "mean"
            for sim_idx = 1:numel(p.outlierparam.meanlist)
                p.outlierparam.mean = p.outlierparam.meanlist(sim_idx);
                output = compute_gnss_ecef(p,eph,obs);
                output.sim_idx = sim_idx;
                outputcell{sim_idx} = output;
            end
        case "count"
            for sim_idx = 1:numel(p.outlierparam.countlist)
                p.outlierparam.count = p.outlierparam.countlist(sim_idx);
                output = compute_gnss_ecef(p,eph,obs);
                output.sim_idx = sim_idx;
                outputcell{sim_idx} = output;
            end
        case "width"
            for sim_idx = 1:numel(p.outlierparam.widthlist)
                p.outlierparam.width = p.outlierparam.widthlist(sim_idx);
                output = compute_gnss_ecef(p,eph,obs);
                output.sim_idx = sim_idx;
                outputcell{sim_idx} = output;
            end
    end
end
end

