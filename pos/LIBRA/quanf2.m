function quan=quanf2(alfa)
n = 30;
p = 4;
rk = p +1;
quan=floor(2*floor((n+rk+1)/2)-n+2*(n-floor((n+rk+1)/2))*alfa);