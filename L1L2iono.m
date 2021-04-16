F = p.L2freq^2/(p.L1freq^2-p.L2freq^2);

L1_sv3 = p.code_bia.GPS.bia_C1C(3);
L2_sv3 = p.code_bia.GPS.bia_C2L(3);
% L1_sv3 = eph.GPS.TGD(3,1);
% L2_sv3 = (p.L1freq^2/p.L2freq^2)*eph.GPS.TGD(3,1);

L1_sv4 = p.code_bia.GPS.bia_C1C(4);
L2_sv4 = p.code_bia.GPS.bia_C2L(4);
% L1_sv4 = eph.GPS.TGD(4,1);
% L2_sv4 = (p.L1freq^2/p.L2freq^2)*eph.GPS.TGD(4,1);
tec_3 = obs.GPS(2).data.P(3,:)-p.c*L2_sv3*1e-9...
    -obs.GPS(1).data.P(3,:)+p.c*L1_sv3*1e-9;
tec_4 = obs.GPS(2).data.P(4,:)-p.c*L2_sv4*1e-9...
    -obs.GPS(1).data.P(4,:)+p.c*L1_sv4*1e-9;
t = p.t(4201:13000);
avgtec = NaN(1,13000-4200);
for i = 4201:13000-400
    avgtec(i-4200) = F * (mean(tec_3(i-200:i+400)) - mean(tec_4(i-200:i+400)));
end
figure
plot(t,avgtec,'.')
hold on
plot(tobs,ustec(3,:)-ustec(4,:),'.')
hold on
plot(tobs,centec(3,:)-centec(4,:),'.')
legend('tec','ustec','cne')
grid on
%%
L1 = eph.GPS.TGD(10,1)+0.3827*1e-9;
L2 = (p.L1freq^2/p.L2freq^2)*eph.GPS.TGD(10,1)-1.4821*1e-9;
tec = obs.GPS(2).data.P(10,:)-p.c*L2...
    -obs.GPS(1).data.P(10,:)+p.c*L1;
avg_tec = NaN(1,6600);
t = p.t(401:7000);
for i = 401:7000
    avg_tec(i-400) = F * mean(tec(i-200:i+400));
end
tobs = datetime(obs_t.tr_prime');
figure
plot(t,avg_tec,'.')
hold on
plot(tobs,ustec(10,:),'.')
hold on
plot(tobs,centec(10,:),'.')
legend('tec','ustec','cne')
grid on

