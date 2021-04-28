function [pos,msr_res,GDOP,nsv,dnsv,nprior,dnprior,postxhat,postxhatcovar,...
    delta_x,Jydiag,by,H_pos,outliervec,outlierbin] = ...
    userposlinear_PVA(p,cpt,grdpos,x_prior,P_prior)
% This is solver for computing user position, receiver clock
% bias, satellite system offset.
% Linear mode
% Measurement selection applied
% Input: 
%       s_pos_ecef: 3-by-N Satellite position in ECEF frame.
%       x0 : 3-by-1 initial interative coordinates in ECEF frame.
%       y: m-by-1 Corrected pseudorange.

%-------------------%

num = length(cpt.corr_range); % number of measurement
if p.eb_outlier == 1
    if p.genOutlier == 1        
        % generate outliers for this simulation
        [cpt.outliervec,cpt.outlierbin] = genoutlier(p.outlierparam,num);
    else
        cpt.outliervec = p.outliervec{p.i}; % load outliers for this dataset
        cpt.outlierbin = p.outlierbin{p.i}; % load outliers for this dataset
    end
else
    cpt.outliervec = zeros(num,1);
    cpt.outlierbin = zeros(num,1);
end

outliervec = cpt.outliervec;
outlierbin = cpt.outlierbin;

% % corrupt measurements with outliers
% cpt.corr_range = cpt.corr_range + cpt.outliervec;
% cpt.dp_range   = cpt.dp_range + cpt.outliervec;

%%%%%%%%%%%%%%%%%%%%%%%% Measurement Selection Algorithms %%%%%%%%%%%%%%%%%
[pos.LS{1},msr_res.LS{1},GDOP.LS,nsv.LS,dnsv.LS,nprior.LS,dnprior.LS,postxhat.LS,postxhatcovar.LS,...
    delta_x.LS,Jydiag.LS,by.LS{1},H_pos.LS{1}]= LSlinear_PVA(p,cpt,grdpos,x_prior.LS,P_prior.LS);

if p.eb_LTS == 1
for idx = 1:p.LTSn
    [pos.LTS{idx},msr_res.LTS{idx},GDOP.LTS(idx),nsv.LTS(idx),dnsv.LTS(idx),nprior.LTS(idx),...
    dnprior.LTS(idx),postxhat.LTS{idx},postxhatcovar.LTS{idx},delta_x.LTS{idx},Jydiag.LTS{idx},by.LTS{idx},H_pos.LTS{idx}]= ...
    LTSlinear_PVA_useallprior(p,cpt,grdpos,x_prior.LTS{idx},P_prior.LTS{idx},p.LTSOption(idx));
end
end

if p.eb_RAPS == 1
for idx = 1:p.RAPSn
    [pos.RAPS{idx},msr_res.RAPS{idx},GDOP.RAPS(idx),nsv.RAPS(idx),dnsv.RAPS(idx),nprior.RAPS(idx),...
    dnprior.RAPS(idx),postxhat.RAPS{idx},postxhatcovar.RAPS{idx},delta_x.RAPS{idx},Jydiag.RAPS{idx},by.RAPS{idx},H_pos.RAPS{idx}]= ...
    RAPSlinear_PVA_useallprior(p,cpt,grdpos,x_prior.RAPS{idx},P_prior.RAPS{idx},p.RAPSEps(idx,:));
end
end

if p.eb_TD == 1
for idx = 1:p.TDn
    [pos.TD{idx},msr_res.TD{idx},GDOP.TD(idx),nsv.TD(idx),dnsv.TD(idx),nprior.TD(idx),...
        dnprior.TD(idx),postxhat.TD{idx},postxhatcovar.TD{idx},delta_x.TD{idx},Jydiag.TD{idx},by.TD{idx},H_pos.TD{idx}]= ...
    TDlinear_PVA_useallprior(p,cpt,grdpos,x_prior.TD{idx},P_prior.TD{idx},p.TDLambda(idx));
end
end

if p.eb_MShb == 1    
for idx = 1:p.MShbn
    [pos.MShb{idx},msr_res.MShb{idx},GDOP.MShb(idx),nsv.MShb(idx),dnsv.MShb(idx),nprior.MShb(idx),...
    dnprior.MShb(idx),postxhat.MShb{idx},postxhatcovar.MShb{idx},delta_x.MShb{idx},Jydiag.MShb{idx},by.MShb{idx},H_pos.MShb{idx}]=...
    MSHuberlinear_PVA_useallprior(p,cpt,grdpos,x_prior.MShb{idx},P_prior.MShb{idx},p.MShbConst(idx));
end
end

if p.eb_MStk == 1
for idx = 1:p.MStkn
    [pos.MStk{idx},msr_res.MStk{idx},GDOP.MStk(idx),nsv.MStk(idx),dnsv.MStk(idx),nprior.MStk(idx),...
    dnprior.MStk(idx),postxhat.MStk{idx},postxhatcovar.MStk{idx},delta_x.MStk{idx},Jydiag.MStk{idx},by.MStk{idx},H_pos.MStk{idx}]= ...
    MSTukeylinear_PVA_useallprior(p,cpt,grdpos,x_prior.MStk{idx},P_prior.MStk{idx},p.MStkConst(idx));
end
end
        
% if p.eb_LSS == 1
% for idx = 1:p.LSSn
%     [pos.LSS{idx},msr_res.LSS{idx},GDOP.LSS(idx),nsv.LSS(idx),dnsv.LSS(idx),nprior.LSS(idx),...
%     dnprior.LSS(idx),postxhat.LSS{idx},postxhatcovar.LSS{idx},delta_x.LSS{idx},Jydiag.LSS{idx},by.LSS{idx},H_pos.LSS{idx}]= ...
%     LSSlinear_PVA(p,cpt,grdpos,x_prior.LSS{idx},P_prior.LSS{idx},p.LSSLambda(idx));
% end
% end


end
