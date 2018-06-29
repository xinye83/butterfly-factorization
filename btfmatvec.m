function [y,nflops] = btfmatvec(B,m,x)
%  [y,nflops] = btfmatvec(B,m,x)
%   mat-vec of a butterfly compressed matrix B and vector x
%
%  Inputs
%     B: struct item --- output from mat2btf
%     m: vector      --- bottom level block sizes
%                        (should be the same as the input m of mat2btf)
%     x:             --- vector where the matrix is applied
%
%  Outputs
%     y:             --- mat-vec
%     nflops:        --- flops count
%
%  Xin Ye, Mar 2017

nb  = length(m);
l   = log2(nb);
nrhs = size(x,2);

m0 = zeros(nb+1,1);
for i=1:nb
    m0(i+1) = m0(i) + m(i);
end

nflops = nb;

xt = cell(nb,1);
yt = xt;

for i=1:nb
    xt{i} = x(m0(i)+1:m0(i+1),:);
    
    xt{i} = xt{i}(B.P{1,1}{1,i}.p,:);
    
    sz1 = size(B.P{1,1}{1,i}.E,1);
    sz2 = size(B.P{1,1}{1,i}.E,2);
    
    xt{i} = xt{i}(1:sz1,:)+B.P{1,1}{1,i}.E*xt{i}(sz1+1:end,:);
    
    nflops = nflops+2*sz1*sz2*nrhs;
end

if l==0
    yt = xt;
else
    for k=1:l
        % at level k
        % 2^(k-1) diagonal blocks, each diagonal is of form [ D_1 ; D_2 ]
        % and D_i is block diagonal with 2^(l-k) blocks
        
        %fprintf('level %d\n',k)
        
        dx = 2^(l-k+1);
        
        for i=1:2^(k-1)
            
            for j=1:dx/2
                yt{j+(i-1)*dx} = [ xt{2*j-1+(i-1)*dx} ; xt{2*j+(i-1)*dx} ];
                yt{j+dx/2+(i-1)*dx} = yt{j+(i-1)*dx};
                
                %size(yt{j+(i-1)*dx})
                %size(B.P{1,k+1}{2*i-1,j}.p)
                
                yt{j+(i-1)*dx} = yt{j+(i-1)*dx}(B.P{1,k+1}{2*i-1,j}.p,:);
                yt{j+dx/2+(i-1)*dx} = yt{j+dx/2+(i-1)*dx}(B.P{1,k+1}{2*i,j}.p,:);
                
                sz1 = size(B.P{1,k+1}{2*i-1,j}.E,1);
                sz2 = size(B.P{1,k+1}{2*i-1,j}.E,2);
                
                yt{j+(i-1)*dx} = yt{j+(i-1)*dx}(1:sz1,:)+...
                    B.P{1,k+1}{2*i-1,j}.E*yt{j+(i-1)*dx}(sz1+1:end,:);
                
                nflops=nflops+2*sz1*sz2*nrhs;
                
                sz1 = size(B.P{1,k+1}{2*i,j}.E,1);
                sz2 = size(B.P{1,k+1}{2*i,j}.E,2);
                
                yt{j+dx/2+(i-1)*dx} = yt{j+dx/2+(i-1)*dx}(1:sz1,:)+...
                    B.P{1,k+1}{2*i,j}.E*yt{j+dx/2+(i-1)*dx}(sz1+1:end,:);
                
                nflops=nflops+2*sz1*sz2*nrhs;
            end
        end
        
        if k==l
            break
        end
        
        xt=yt;
    end    
end

y = zeros(size(x));
for i=1:nb
    sz1 = size(B.U{i,1},1);
    sz2 = size(B.U{i,1},2);
    
    y(1+m0(i):m0(i+1),:) = B.U{i,1}*yt{i};
    
    nflops = nflops+sz1*(2*sz2-1)*nrhs;
end