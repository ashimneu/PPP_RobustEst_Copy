function output = func_Gamma_v(T,n)
In = eye(n);
output = blkdiag((T^(5/2)/sqrt(20)).*In,...
         (T^(3/2)/sqrt(3)).*In, sqrt(T).*In);
end