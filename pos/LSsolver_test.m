function [xk,res_R] = LSsolver_test(p,xp,xv,cpt)
 
y_R = cpt.corr_range;
y_D = cpt.dp_range;
num = length(y_R); % The number of measurement
num_clk = sum(cpt.num_sv~=0);
H = zeros(num,3);
H_D = zeros(num,3);
H_clk = zeros(num,num_clk);
ind = find(cpt.num_sv~=0);
start = 1;
for i = 1:num_clk
    H_clk(start:start+cpt.num_sv(ind(i))-1,i)=1;
    start = start + cpt.num_sv(ind(i));
end
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
for iter=1:p.Nls
    for j=1:num
        R(j)=norm(s_pos_ecef(:,j)-xp(1:3));
        H(j,:) = (xp(1:3)-s_pos_ecef(:,j))'/R(j)+...
            [-s_pos_ecef(2,j)*p.omge/p.c s_pos_ecef(1,j)*p.omge/p.c 0];
        r(j) = R(j)+sagnac(p,s_pos_ecef(:,j),xp);
        ind = find(H_clk(j,:)==1);
        clk_bia(j) = xp(3+ind);
    end
    H_os = [H,H_clk];
    res_R = y_R - r - clk_bia;
    delta_x = (H_os'*H_os)^(-1)*H_os'*res_R;
    xp=xp+delta_x;
    if (norm(delta_x) < p.LSthrsh)
        break;
    end
end
for iter=1:p.Nls
    for j=1:num
        rv(j) = -H(j,:)*(s_v_ecef(:,j)-xv(1:3))+...
            sagnac_v(p,[s_pos_ecef(:,j);s_v_ecef(:,j)],[xp(1:3);xv(1:3)]);
    end
    H_os = [H,ones(num,1)];
    res_D = (y_D - (rv+xv(end)));
    delta_x = (H_os'*H_os)^(-1)*H_os'*res_D;
    xv=xv+delta_x;
    if (norm(delta_x) < p.LSthrsh)
        break;
    end
end
xk = [xp(1:3);xv(1:3);xp(4:end);xv(4)];
end