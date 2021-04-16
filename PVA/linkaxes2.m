function linkaxes2(option)
% INPUT - option.ax
%       - option.axes2link

    if isfield(option,'axes2link')
        if ~strcmp(option.axes2link,'') | ~isempty(option.axes2link)
            ax = option.ax(isvalid(option.ax)); % axes are invalid  if corresponding figure windows are closed.
            linkaxes(option.ax,option.axes2link)
        end
    end
end

