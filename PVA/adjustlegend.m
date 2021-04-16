function [] = adjustlegend(plt,location)
% plt - plot
% set legend orientation to horizontal
% Legend position is set to top left outside of the axis box.
% location of legend compared to axes.

leg = plt.Legend;
set(leg,'Orientation','horizontal');
pltPos = plt.Position;
newPos = leg.Position;

% default - top left & outside of axes.
newPos(1) = pltPos(1); newPos(2) = pltPos(2) + pltPos(4);
set(leg,'Position',newPos);

if nargin == 2
    switch lower(location)
        case 'bottomleft'
        % bottom left & outside of axes.
        newPos(1) = pltPos(1); newPos(2) = pltPos(2) - 0.2*pltPos(4);
        set(leg,'Position',newPos);
    end
end
    
end

