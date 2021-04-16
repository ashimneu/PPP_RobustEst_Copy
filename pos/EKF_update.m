function [p,res_R] = EKF_update(p,cpt)
for kk = 1
ind_s = find(cpt.num_sv([1,3,4]) ~= 0);
xk = [p.state_PVA(1:9);p.state_PVA(9+ind_s);p.state_PVA(13)];
% Pcov = [p.P_cov(1:9,1:9),p.P_cov(1:9,[9+ind_s,13]);
%         p.P_cov([9+ind_s,13],1:9),p.P_cov([9+ind_s,13],[9+ind_s,13])];
ind_no = find(cpt.num_sv([1,3,4]) == 0);
Pcov = p.P_cov;
Pcov(9+ind_no,:) = [];
Pcov(:,9+ind_no) = [];
y_R = cpt.corr_range;
y_D = cpt.dp_range;
num = length(y_R); % The number of measurement
num_clk = sum(cpt.num_sv~=0);
H = zeros(num,3);
H_clk = zeros(num,num_clk);
ind = find(cpt.num_sv~=0);
start = 1;
for i = 1:num_clk
    H_clk(start:start+cpt.num_sv(ind(i))-1,i)=1;
    start = start + cpt.num_sv(ind(i));
end
H_clk(:,1) = 1;
R = zeros(num,1);  
r = zeros(num,1);
rv = zeros(num,1);
clk_bia = zeros(num,1);
lamda = p.c*[ones(cpt.num_sv(1),1)/p.L1freq;ones(cpt.num_sv(2),1)/p.L1freq;...
    ones(cpt.num_sv(3),1)/p.E1freq;ones(cpt.num_sv(4),1)/p.B1freq];
y_D = -lamda.*y_D;
if p.post_mode == 1
    if p.IGS_enable == 1
        s_pos_ecef = cpt.s_pos_prc;
        s_v_ecef = cpt.s_v_ecef;
    else
        s_pos_ecef = cpt.s_pos_ecef;
        s_v_ecef = cpt.s_v_ecef;
    end
else
    s_pos_ecef = cpt.s_pos_ecef;
    s_v_ecef = cpt.s_v_ecef;
end
% if p.mk == 1
%     xk = [xk;0]; %Add Iono factor term for PPP
% end
dx_list = nan(size(xk,1),5); xk_list = nan(size(xk,1),5); 
res_list = nan(2*num,5); msr_res_list = nan(2*num,5);
if kk == 1
    iter_count = 1; 
elseif kk == 2
    iter_count = 1; 
end
    
for iter=1:iter_count%5
for j=1:num
    R(j)=norm(s_pos_ecef(:,j)-xk(1:3));
    H(j,:) = (xk(1:3)-s_pos_ecef(:,j))'/R(j)+...
        [-s_pos_ecef(2,j)*p.omge/p.c s_pos_ecef(1,j)*p.omge/p.c 0];
    rv(j) = -H(j,:)*(s_v_ecef(:,j)-xk(4:6))+...
        sagnac_v(p,[s_pos_ecef(:,j);s_v_ecef(:,j)],xk);
    r(j) = R(j)+sagnac(p,s_pos_ecef(:,j),xk);
    ind = find(H_clk(j,:)==1);
    clk_bia(j) = sum(xk(9+ind));
end
H_os = [H,zeros(size(H)),zeros(num,3),H_clk,zeros(num,1);
    zeros(size(H)),H,zeros(num,3),zeros(size(H_clk)),ones(num,1)];
res_R = y_R - r - clk_bia;
res_D = (y_D - (rv+xk(end)));
%  Innovation Covariance
S = H_os*Pcov*H_os'+eye(num*2);
% Optimal Kalman Gain
Kk = Pcov*H_os'*S^(-1);
% Posterior state estimate
dx = Kk*[res_R;res_D];
xk = xk + dx;
msr_res = [res_R; res_D] - H_os*dx;
dx_list(:,iter) = dx; xk_list(:,iter) = xk; 
res_list(:,iter) = [res_R; res_D]; msr_res_list(:,iter) = msr_res;
p.msr_res = msr_res(1:num);
% Posterior state covariance estimate
T = Kk*H_os;
Pcov = (eye(size(T))-T)*Pcov;
% end

p.state_PVA([1:9,9+ind_s,13]) = xk;
ivec = [1:9,9+ind_s,13];
for l = 1:length(ivec)
    p.P_cov(ivec(l),[1:9,9+ind_s,13]) = Pcov(l,:);
end

% if kk == 2 %~(p.i == 26)
%     p.state_PVA([1:9,9+ind_s,13]) = xk;
%     ivec = [1:9,9+ind_s,13];
%     for l = 1:length(ivec)
%         p.P_cov(ivec(l),[1:9,9+ind_s,13]) = Pcov(l,:);
%     end
% else
% %     p.state_PVA([9+ind_s,end]) = [xk(9+ind_s); xk(end)];
%     p.state_PVA(10) = xk(10);
% end
% p.P_cov = Pcov;
end % iter

% if kk == 1
%     dv = dx_list(4:6,:); 
%     v = xk_list(4:6,:);
%     dv15 = xk_list(4:6,1) - xk_list(4:6,end);
%     p.vnorm = norm(dv15);
% end

end %for kk


end