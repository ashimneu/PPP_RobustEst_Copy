function [newA,idx_nonzerocols,idx_nonzerorows] = remove0colrow(A,option)
% Removes rows/columns of A that consist of only 0 entires.

idx_nonzerorows = nan(size(A,2));
idx_nonzerocols = nan(size(A,1));

switch lower(option)
    case "column"
        newA = A;
        zerocol = ~any(A,1);    % assigns 0 if col of 0s is present
        newA(:,zerocol) = [];   % removes cols of 0s
        idx_nonzerocols = find(~zerocol);
        
    case "row"
        newA = A;
        zerorow = ~any(newA,2); % assigns 0 if row of 0s is present
        newA(zerorow,:) = [];   % removes rows of 0s
        idx_nonzerorows = find(~zerorow);
        
    case "both"
        newA = A;
        zerocol = ~any(A,1);    % assigns 0 if col of 0s is present
        newA(:,zerocol) = [];   % removes cols of 0s
        zerorow = ~any(newA,2); % assigns 0 if row of 0s is present
        newA(zerorow,:) = [];   % removes rows of 0s
        idx_nonzerocols = find(~zerocol);
        idx_nonzerorows = find(~zerorow);
end
end

