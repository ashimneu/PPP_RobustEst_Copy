function [x_out,P_out] = adjust_entries(p,x_in,P_in,option)
% function [x_out,P_out,Pbx_out] = adjust_entries(p,x_in,P_in,option,Pbx_in)

x_out   = x_in;
P_out   = P_in;

% active sat constellation besides GPS at current measurement epoch
ind = p.ind; 
      
switch option
    case 0 
        % do nothing, make no adjustments
    case 1
        if ~isempty(ind)
            ISB = [p.ISBglo;p.ISBgal;p.ISBbds];
            ISB_cov = [p.ISBglo_cov;p.ISBgal_cov;p.ISBbds_cov];
            ISB_active = ISB(ind);
            ISB_cov_active = ISB_cov(ind);

            x_out = [x_in(1:10); ISB_active; x_in(10+1:end)];
            P_tmp = blkdiag(P_in(1:10,1:10), diag(ISB_cov_active), P_in(11,11));
            P_tmp(1:10,end) = P_in(1:10,end);
            P_tmp(end,1:10) = P_in(end,1:10);
            P_out = P_tmp;
        end
    case 2
        if ~isempty(ind)
            x_out = [x_in(1:10); x_in(end)];
            P_tmp = blkdiag(P_in(1:10,1:10), P_in(end,end));
            P_tmp(1:10,end) = P_in(1:10,end);
            P_tmp(end,1:10) = P_in(end,1:10);
            P_out = P_tmp;
        end
%         if nargout == 3 && nargin == 5
%             Pbx_tmp = blkdiag(Pbx_in(1:10,1:10), Pbx_in(end,end));
%             Pbx_tmp(1:10,end) = Pbx_in(1:10,end);
%             Pbx_tmp(end,1:10) = Pbx_in(end,1:10);
%             Pbx_out = Pbx_tmp;
%         end
end


    
end
