BDS.t_oc = cell(65,1);
tget = str2double(strsplit("2021 02 27 08 00 00"));
tget = datetime(tget)+seconds(14); % seconds diff from GPS time to BDS time
tget = [tget.Year,tget.Month,tget.Day,tget.Hour,tget.Minute,tget.Second];
indx = 1;
prn = 20;
[~,~,BDS.t_oc{prn}(indx)] = date2gpst(tget); % Represent BDS time by GPS time
data = sscanf("-9.207843104377e-04-5.495159882685e-12 0.000000000000e+00",'%f');
BDS.a_f0(prn,indx) = data(1); % SV clock bias (seconds)
BDS.a_f1(prn,indx) = data(2); % SV clock drift (sec/sec)
BDS.a_f2(prn,indx) = data(3); % SV clock drift rate (sec/sec2)

data = sscanf("1.000000000000e+00-1.731875000000e+02 3.720869274855e-09 5.150960860367e-01",'%f');
BDS.AODE(prn,indx) = data(1); % Age of Data, Ephemeris
BDS.C_rs(prn,indx) = data(2); % Crs (meters)
BDS.Delta_n(prn,indx) = data(3); % Delta n (radians/sec)
BDS.M_0(prn,indx) = data(4); % M0 (radians)

data = sscanf("-8.486211299896e-06 7.530053844675e-04 9.120907634497e-06 5.282637945175e+03",'%f');
BDS.C_uc(prn,indx) = data(1); % Cuc (radians)
BDS.e(prn,indx) = data(2); % e Eccentricity
BDS.C_us(prn,indx) = data(3); % Cus (radians)
BDS.sqrtA(prn,indx) = data(4); % sqrt(A) (sqrt(m))

data = sscanf("5.472000000000e+05 4.284083843231e-08-4.089488597329e-01-1.536682248116e-08",'%f');
BDS.t_oe(prn,indx) = data(1); % Toe Time of Ephemeris (sec of BDS week)
BDS.IODE(prn,indx) = mod(floor(data(1)/720),240);
BDS.C_ic(prn,indx) = data(2); % Cic (radians)
BDS.Omega_0(prn,indx) = data(3); % OMEGA0 (radians)
BDS.C_is(prn,indx) = data(4); % Cis (radians)

data = sscanf("9.657224475648e-01 1.812343750000e+02-7.765824654454e-01-6.890287008111e-09",'%f');
BDS.i_0(prn,indx) =  data(1); % i0 (radians)
BDS.C_rc(prn,indx) = data(2); % Crc (meters)
BDS.Omega(prn,indx) = data(3); % omega (radians)
BDS.OmegaDot(prn,indx) = data(4); % OMEGA DOT (radians/sec)

data = sscanf("-4.603763193788e-10 0.000000000000e+00 7.900000000000e+02 0.000000000000e+00",'%f');
BDS.IDOT(prn,indx) = data(1); % IDOT (radians/sec)
BDS.week_num(prn,indx) = data(3); % BDT Week #

data = sscanf("2.000000000000e+00 0.000000000000e+00 2.300000000000e-08 2.300000000000e-08",'%f');
BDS.SV_acc(prn,indx) = data(1); % SV accuracy (meters) See BDS ICD 200H Section 20.3.3.3.1.3
BDS.SV_health(prn,indx) = data(2); % SV health
BDS.TGD1(prn,indx) = data(3); % TGD1   B1/B3          (seconds)
BDS.TGD2(prn,indx) = data(4); % TGD2   B2/B3          (seconds)

data = sscanf("5.472000000000e+05 1.000000000000e+00",'%f');
BDS.trans_time(prn,indx) = data(1); % Transmission time of message
BDS.ADOC(prn,indx) = data(2); % Age of Data Clock

%%
tget = str2double(strsplit("2021 02 27 10 00 00"));
GPS.t_oc = cell(65,1);
[~,~,GPS.t_oc{prn}(indx)] = date2gpst(tget); % Time of broadcast (seconds)
data = sscanf("5.244966596365e-04-2.273736754432e-13 0.000000000000e+00",'%f');
GPS.a_f0(prn,indx) = data(1); % SV clock bias (seconds)
GPS.a_f1(prn,indx) = data(2); % SV clock drift (sec/sec)
GPS.a_f2(prn,indx) = data(3); % SV clock drift rate (sec/sec2)

data = sscanf("7.500000000000e+01-1.058437500000e+02 5.474513749758e-09-6.397341492094e-01",'%f');
GPS.IODE(prn,indx) = data(1); % Issue of Data, Ephemeris
GPS.C_rs(prn,indx) = data(2); % Crs (meters)
GPS.Delta_n(prn,indx) = data(3); % Delta n (radians/sec)
GPS.M_0(prn,indx) = data(4); % M0 (radians)

data = sscanf("-5.505979061127e-06 5.841623409651e-03 2.807006239891e-06 5.153675813675e+03",'%f');
GPS.C_uc(prn,indx) =  data(1); % Cuc (radians)
GPS.e(prn,indx) = data(2); % e Eccentricity
GPS.C_us(prn,indx) = data(3); % Cus (radians)
GPS.sqrtA(prn,indx) = data(4); % sqrt(A) (sqrt(m))

data = sscanf("5.544000000000e+05-3.725290298462e-09-9.159405288281e-01-2.048909664154e-08",'%f');
GPS.t_oe(prn,indx) = data(1); % Toe Time of Ephemeris (sec of GPS week)
GPS.C_ic(prn,indx) = data(2); % Cic (radians)
GPS.Omega_0(prn,indx) = data(3); % OMEGA0 (radians)
GPS.C_is(prn,indx) = data(4); % Cis (radians)

data = sscanf("9.370117791017e-01 3.112500000000e+02 2.866666821757e+00-8.894656212419e-09",'%f');
GPS.i_0(prn,indx) =  data(1); % i0 (radians)
GPS.C_rc(prn,indx) = data(2); % Crc (meters)
GPS.Omega(prn,indx) = data(3); % omega (radians)
GPS.OmegaDot(prn,indx) = data(4); % OMEGA DOT (radians/sec)

data = sscanf("-2.792973481413e-10 1.000000000000e+00 2.146000000000e+03 0.000000000000e+00",'%f');
GPS.IDOT(prn,indx) = data(1); % IDOT (radians/sec)
GPS.CodesOnL2(prn,indx) = data(2); % Codes on L2 channel
GPS.week_num(prn,indx) = data(3); % GPS Week #
GPS.L2Pflag(prn,indx) = data(4); % L2 P data flag

data = sscanf("2.000000000000e+00 0.000000000000e+00-8.381903171539e-09 7.500000000000e+01",'%f');
GPS.SV_acc(prn,indx) = data(1); % SV accuracy (meters) See GPS ICD 200H Section 20.3.3.3.1.3
GPS.SV_health(prn,indx) = data(2); % SV health
GPS.TGD(prn,indx) = data(3); % TGD (seconds)
GPS.IODC(prn,indx) = data(4); % IODC Issue of Data, Clock

GPS.trans_time(prn,indx) = 0; % Transmission time of message
GPS.fit_interval(prn,indx) = 0; % Fit Interval in hours

%%
obs = [];
tidx = 1;
[rt.week, rt.dow, rt.sow] = date2gnsst([2021 2 27 8 37 18]);
rt.sow = round(rt.sow)- p.bds.lps_gps;
sys = 'BDS';
p.post_mode = 0;
[tp,dt_sv,sat] = inverse_eph(p,BDS,obs,prn,tidx,rt.sow,sys,p.P_base);
r = norm(p.P_base-sat.pos_ecef)+sagnac(p,sat.pos_ecef,p.P_base)-dt_sv*p.c

