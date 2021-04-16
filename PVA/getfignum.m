function figno = getfignum(figtag)
%   figtag - user specified figure tag
    destFig = findobj('type','figure','tag',figtag);
    if ~(isempty(destFig))
        figno = destFig.Number; % figure number
    else
        fig =  figure('tag',figtag); 
        figno = fig.Number;
    end
end