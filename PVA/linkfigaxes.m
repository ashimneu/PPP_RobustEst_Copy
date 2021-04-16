function [] = linkfigaxes(fignum,linkaxis)
ax = [];
for ii = 1:numel(fignum)
    if ishandle(ii)
        fig = figure(ii);
        ax = [ax; findall(fig,'type','axes')];
    end
end
linkaxes(ax,linkaxis);
end

