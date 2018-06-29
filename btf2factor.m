function [U,P] = btf2factor(B)
% convert the struct item B into factors U,P{1},P{2},...,P{l+1}
%     A = U * P{l+1} * P{l} * ... * P{1}
% (used for debugging)

U = [];
for i=1:length(B.U)
    U = blkdiag(U,B.U{i,1});
end

P = cell(length(B.P),1);

P{1} = [];
for i=1:length(B.P{1,1})
    P0 = B.P{1,1}{1,i}.E;
    P0 = [eye(size(P0,1)),P0];
    P0 = P0(:,invp(B.P{1,1}{1,i}.p));
    
    P{1} = blkdiag(P{1},P0);
end

% figure(1)
% spy(P{1})
% title('P{1}')

l = length(B.P)-1;

for k=1:l
    P{k+1} = [];
    
    for i=1:2^(k-1)
        D1 = [];
        D2 = [];
        
        for j=1:2^(l-k)
            P0 = B.P{1,k+1}{2*i-1,j}.E;
            P0 = [eye(size(P0,1)),P0];
            P0 = P0(:,invp(B.P{1,k+1}{2*i-1,j}.p));
            
            D1 = blkdiag(D1,P0);
            
            P0 = B.P{1,k+1}{2*i,j}.E;
            P0 = [eye(size(P0,1)),P0];
            P0 = P0(:,invp(B.P{1,k+1}{2*i,j}.p));
            
            D2 = blkdiag(D2,P0);
        end
        
        P{k+1} = blkdiag(P{k+1},[D1;D2]);
    end
    
%     figure(k+1)
%     spy(P{k+1})
%     title(['P{',int2str(k+1),'}'])
end

% figure(l+2)
% spy(U)
% title('U')


