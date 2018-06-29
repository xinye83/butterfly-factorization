clear all
close all

e = 6;


N = 100*2^e;
n = N/2;
fprintf('matrix size:\t\t%d\n\n',n)

A = dft((1:n),(n+1:N),N);
x = randn(n,1);
y = A*x;

tol = 1e-8;
m   = 200*ones(n/200,1);

[B,nflops1] = mat2btf(A,m,tol);

% computing residual for the matrix
% only use this for small problems
res = factorres(A,B);
fprintf('compression res:\t%e\n',res)

fprintf('compression cost:\t%e\n',nflops1)

fact=whos('B');
fprintf('memory:\t\t\t%f MB\n\n',fact.bytes/1024/1024)

[y1,nflops2] = btfmatvec(B,m,x);

fprintf('multiplication res:\t%e\n',norm(y-y1)/norm(y))
fprintf('multiplication cost:\t%e\n',nflops2)