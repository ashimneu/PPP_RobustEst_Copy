function output = func_F_v(T,lambda_a,n)
a1 = exp(-lambda_a*T);
a2 = (1 - a1)/lambda_a;
a3 = (lambda_a*T - 1 + a1)/(lambda_a^2);
In = eye(n);
Zn = zeros(n,n);
output = [In, T*In, a3*In;...
          Zn, In,   a2*In;...
          Zn, Zn,   a1*In];
end