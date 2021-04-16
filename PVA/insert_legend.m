function output_array = insert_legend(cell_val, cell_array, idx)
    if idx == -1
        output_array = cat(1,cell_array,cell_val); % add to end of array
    elseif size(cell_array,1) == 1
        output_array = cat(1,cell_val,cell_array); % insert at index "idx" of array
    else
        output_array = cat(1,cell_array(1:idx-1),cell_val,cell_array(idx:end));
    end
end