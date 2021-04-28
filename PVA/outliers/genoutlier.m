function [outlier_vector,outlier_binary] = genoutlier(outlierparam,num)

    if (outlierparam.mean -  outlierparam.width) < 0 
        lbound = 0;
    else
        lbound =  outlierparam.mean -  outlierparam.width;
    end
    ubound =  outlierparam.mean +  outlierparam.width;
    outlier_vector = zeros(num,1);
    outlier_vector(1:outlierparam.count) = lbound + (ubound-lbound)*rand(outlierparam.count,1); 
    outlier_vector = outlier_vector(randperm(num)); % permutates entries in the vector
    outlier_binary = any(outlier_vector,2); % true binary vector that reflects presence of outlier in the vector
end

