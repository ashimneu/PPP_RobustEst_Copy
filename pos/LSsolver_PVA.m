function [xk,res_R] = LSsolver_PVA(p,xk,cpt)
 
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
% for iter=1:p.Nls
    for j=1:num
        R(j)=norm(s_pos_ecef(:,j)-xk(1:3));
        H(j,:) = (xk(1:3)-s_pos_ecef(:,j))'/R(j)+...
            [-s_pos_ecef(2,j)*p.omge/p.c s_pos_ecef(1,j)*p.omge/p.c 0];
        rv(j) = -H(j,:)*(s_v_ecef(:,j)-xk(4:6))+...
                    sagnac_v(p,[s_pos_ecef(:,j);s_v_ecef(:,j)],xk);
        r(j) = R(j)+sagnac(p,s_pos_ecef(:,j),xk);
        ind = find(H_clk(j,:)==1);
        clk_bia(j) = sum(xk(6+ind));
    end
    H_os = [H,zeros(size(H)),H_clk,zeros(num,1);
            zeros(size(H)),H,zeros(size(H_clk)),ones(num,1)];
    res_R = y_R - r - clk_bia;
    res_D = (y_D - (rv+xk(end)));
delta_x = (H_os'*H_os)^(-1)*H_os'*[res_R;res_D];
xk=xk+delta_x;
% if (norm(delta_x([1:3,7])) < p.LSthrsh)
%      break;
% end 
% if (iter>p.Nls)&& (norm(delta_x) > p.LSthrsh)
%      warning('Postion path length iteration failed in user_pos calculation');
% end  

% end
% GDOP = sqrt(trace((H_os'*H_os)^(-1)));
% if GDOP_check(p,GDOP)==1  
%     pos = xk(1:3);
%     clock_bias = xk(4);
% else
%     pos = [];
%     clock_bias = [];
% end      

end