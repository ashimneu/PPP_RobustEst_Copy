function output_array = insert_vector(val, array, idx)
    if idx == -1
        output_array = cat(2,array,val); % add to end of array
    elseif size(array,1) == 1
        output_array = cat(2,[],val,array(idx:end)); % insert at index "idx" of array
    else
        output_array = cat(2,array(:,1:idx-1),val,array(:,idx:end));
    end
end