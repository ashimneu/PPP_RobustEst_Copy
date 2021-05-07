function [ned_err_solver,hor_err_solver,err_solver,GDOP_solver,nsv,dnsv,nprior,dnprior,msr_res_solver,by_solver,H_pos_solver,res_std_solver] =...
    save_errNorm_res_GDOP...
    (solver_n,Grdpose,re_pos_solver,GDOP,nsv,dnsv,nprior,dnprior,max_num_sv,ind_mark,msr_res_cells,by_cells,H_pos_cells,res_std_cells)
    
        
    err_solver     = NaN(solver_n,1); 
    hor_err_solver = NaN(solver_n,1);    
    ned_err_solver = NaN(solver_n,1);
    GDOP_solver    = NaN(solver_n,1);

    msr_res_solver = [];
    by_solver      = [];
    H_pos_solver   = [];
    res_std_solver = [];
    
%     msr_res_solver = NaN(max_num_sv*solver_n,1);
%     by_solver      = NaN(max_num_sv*solver_n,1);
%     H_pos_solver   = NaN(max_num_sv*solver_n,3);

    for idx=1:solver_n
        [pos_llh,~,~]=ecef2llh_iter(re_pos_solver{idx});
        R_e2g =ll2R(pos_llh); % rotation matrix from ecef 2 geodetic frame
        err_pos = Grdpose - re_pos_solver{idx}; % position err = true - estimated position
        ned_err=R_e2g*err_pos;    
        ned_err_solver(idx) = norm(ned_err);
        hor_err_solver(idx) = norm(ned_err(1:2));        
        err_solver(idx)     = norm(err_pos);
        GDOP_solver(idx)    = GDOP(idx);
        %------------------------------%
        num = sum(ind_mark);
        tmp_res_R = NaN(max_num_sv,1); tmp_res_D = NaN(max_num_sv,1);
        tmp_res_R(ind_mark) = msr_res_cells{idx}(1:num); tmp_res_D(ind_mark) = msr_res_cells{idx}(num+1:end);
        msr_res_solver = [msr_res_solver; tmp_res_R; tmp_res_D];
        %------------------------------%
        tmp_by_R = NaN(max_num_sv,1); tmp_by_D = NaN(max_num_sv,1);
        tmp_by_R(ind_mark) = by_cells{idx}(1:num); tmp_by_D(ind_mark) = by_cells{idx}(num+1:end);
        by_solver = [by_solver; tmp_by_R ;tmp_by_D];
        %------------------------------%
        tmp_H_pos_R = NaN(max_num_sv,6); tmp_H_pos_D = NaN(max_num_sv,6);
        tmp_H_pos_R(ind_mark,:) = H_pos_cells{idx}(1:num,:); tmp_H_pos_D(ind_mark,:) = H_pos_cells{idx}(num+1:end,:);
        H_pos_solver = [H_pos_solver; tmp_H_pos_R; tmp_H_pos_D];
        %------------------------------%
        tmp_res_std_R = NaN(max_num_sv,1); tmp_res_std_D = NaN(max_num_sv,1);
        tmp_res_std_R(ind_mark) = res_std_cells{idx}(1:num); tmp_res_std_D(ind_mark) = res_std_cells{idx}(num+1:end);
        res_std_solver = [res_std_solver; tmp_res_std_R; tmp_res_std_D];
    end
end

