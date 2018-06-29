function res = factorres(A,B)

[U,P] = btf2factor(B);

A1=U;

for i=1:length(P)
    A1 = A1*P{end-i+1};
end

res = norm(A-A1)/norm(A);