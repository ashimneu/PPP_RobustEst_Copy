function linkaxes2(option)
% INPUT - option.ax
%       - option.axes2link

    if isfield(option,'axes2link')
        if ~strcmp(option.axes2link,'') | ~isempty(option.axes2link)
            linkaxes(option.ax,option.axes2link)
        end
    end
end

