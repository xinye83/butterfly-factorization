function M = dft(J,K,n)
% Generate the submatrix of DFT matrix of order n
%  with row index J and column index K

J = reshape(J,length(J),1);
K = reshape(K,length(K),1);

M = J.*K'; 
% The above line only works for R2016b and later, use the following line
% if you have earlier version of Matlab
% M = bsxfun(@times,J,K');

M = exp(2*pi*1i*M/n);