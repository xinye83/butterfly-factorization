function [B,nflops] = mat2btf(A,m,tol)
%  [B,nflops] = mat2btf(A,m,tol)
%   Compress matrix A into a butterfly factorization [O'Neil Rokhlin 2007]
%
%  Inputs:
%      A: n x n matrix         --- the matrix to be compressed
%      m: vector of length 2^l --- bottom level block sizes
%    tol: real                 --- compression tolerance
%
%  Outputs:
%        B: struct item        --- contains all factors in factorization
%           A = B.U * B.P{1} * B.P{2} * ... * B.P{l+1}
%   nflops: integer            --- flops count
%
%  Xin Ye, Mar 2017

n  = size(A,1); % matrix size
nb = length(m); % # of blocks

if nb~=pow2(nextpow2(nb))
    error('Length of m must be a power of 2.');
end

l  = log2(nb);  % # of levels

m0 = zeros(nb+1,1);
for i=1:nb
    m0(i+1) = m0(i) + m(i);
end

if n~=m0(end)
    error('Sum of m does not equal to n.');
end

B.U = cell(nb,1);
B.P = cell(1,l+1);
nflops = nb;

Ilast = cell(1,nb);

for i=1:nb
    Ilast{1,i} = (m0(i)+1:m0(i+1))';
end

for k=0:l
    % at level k
    % i=1:2^k, j=1:2^(l-k), row size=n/2^k
    B.P{k+1} = cell(2^k,2^(l-k));
    Ithis = cell(2^k,2^(l-k));
    
    di = 2^(l-k); dj = 2^k;
    
    %fprintf('level %d of %d\n',k,l)
    
    for i=1:2^k
        indi = (m0((i-1)*di+1)+1:m0(i*di+1))';
        
        for j=1:2^(l-k)
            
            %fprintf('(%d,%d) block\n',i,j)
            
            if k==0
                indj = Ilast{1,j};
            else
                indj = [ Ilast{ceil(i/2),2*j-1} ;...
                    Ilast{ceil(i/2),2*j} ];
            end
            
            [J,E,p,nflops0] = compr2id(A(indi,indj).',tol);
            
            B.P{k+1}{i,j}.E = E.';
            B.P{k+1}{i,j}.p = p;
            Ithis{i,j} = indj(J); % global indices
            
            nflops = nflops+nflops0;
            
        end
        
    end
    
    if k==l
        for i=1:nb
            B.U{i,1} = A( (m0(i)+1:m0(i+1))' , Ithis{i,1} );
        end
        
        break
    end
    
    Ilast = Ithis;
end