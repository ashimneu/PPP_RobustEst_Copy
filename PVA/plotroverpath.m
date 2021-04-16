Grdpos = p.Grdpos.pos;
p = output.p;
wgs84 = wgs84Ellipsoid(); % meters


N = size(p.Grdpos.pos,2); % total epochs for which true pos are available
Tstart = 1; % 1130;    % [seconds] rover start epoch
Tend   = N; % [seconds] rover end epoch
traverse_epoch = Tstart:1:Tend;
p_grd = Grdpos(:,traverse_epoch); % rover true position
grd0 = ecef2lla(Grdpos(:,1)'); % [degrees] initial position
[trueN,trueE,trueD] = ecef2ned(p_grd(1,:),p_grd(2,:),...
    p_grd(3,:),grd0(1),grd0(2),grd0(3),wgs84);

solvername = "ls";
figtag = strcat("traj2d_",solvername); fignum = getfignum(figtag);
figure(fignum); clf; hold on; grid on
plot(trueE,trueN, 'k-')

xlabel('East, meters')
ylabel('North, meters')

propsteps = output.p.numPropSteps;
propTstart = (Tstart-1) * propsteps + 1; % [seconds]
propTend   = (Tend-1) * propsteps;       % [seconds]
p_prop = output.prior_state.LS(1:3,propTstart:propTend); % rover pos estimates after each propagation step
[prN,prE,prD] = ecef2ned(p_prop(1,:),p_prop(2,:),p_prop(3,:),grd0(1),grd0(2),grd0(3),wgs84);
plot(prE,prN,'g.')

p_post = output.post_state.LS(1:3,traverse_epoch); % rover pos estimates after measurement update
[psN,psE,psD] = ecef2ned(p_post(1,:),p_post(2,:),p_post(3,:),grd0(1),grd0(2),grd0(3),wgs84);
plot(psE,psN, 'm.')

startmark = scatter(trueE(1),trueN(1),'sc','HandleVisibility','off');
scatter(trueE(end),trueN(end),'pc','HandleVisibility','off');
text(trueE(1),trueN(1)-50,'start');
text(trueE(end),trueN(end),'end');

legend('RTK/True','KF Prop.','LS')
xlim([-2300 2300])
ylim([-700 3800])


%%
% solvername = "ls";
figtag = strcat("errxyz2d_",solvername); fignum = getfignum(figtag);
figure(fignum); clf; 
subplot(311); hold on; grid on
plot(traverse_epoch, psE - trueE, '.')
title('East')
ylabel('$\bf{\hat{p}_E} - \bf{p}_E$, meters','Interpreter','latex')
xlabel('GPS time in seconds')

subplot(312); hold on; grid on
plot(traverse_epoch, psN - trueN, '.')
title('North')
ylabel('$\bf{\hat{p}_N} - \bf{p}_N$, meters','Interpreter','latex')
xlabel('GPS time in seconds')

subplot(313); hold on; grid on
plot(traverse_epoch,psD - trueD,'.')
title('Depth')
ylabel('$\bf{\hat{p}_D} - \bf{p}_D$, meters','Interpreter','latex')
xlabel('GPS time in seconds')
sgtitle('Measurement update using Least Squares')