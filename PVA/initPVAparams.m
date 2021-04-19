function p = initPVAparams(p)
    n = p.n;
    p.x_prior = [p.x_v_prior; p.x_c_prior];
    p.P_prior = diag([p.pva_cov_prior; p.clk_cov_prior]);
    
    % PVA State Space model matrices
    F_v   = func_F_v(p.T,p.lam_a,n);
    Gam_v = func_Gamma_v(p.T,n);
    Q_v   = func_Q_v(p.sig_p,p.sig_v,p.sig_a,n);

    F_c   = func_F_c(p.T,p.lam_cdrift);
    Gam_c = func_Gamma_c(p.T);
    Q_c   = func_Q_c(p.sig_cbias,p.sig_ISB_E,p.sig_ISB_B,p.sig_cdrift);
    
    p.F   = blkdiag(F_v,F_c);
    p.Gam = blkdiag(Gam_v,Gam_c);
    p.Q   = blkdiag(Q_v,Q_c);

end