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