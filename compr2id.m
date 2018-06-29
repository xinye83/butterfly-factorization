function [J,E,p,nflops] = compr2id(A,tol)

% [~,R,nflops,r,p] = mgsclpv0(A','tol',tol);
[~,R,nflops,r,p] = mgsclpv(A',tol);

J = p(1:r); I = p(r+1:end);
E = R(:,J)\R(:,I); % in mgsclpv, cols of R are not actully permuted
E = E';

% [U,S] = svd(A,'econ');
% 
% sv = diag(S);
% r = length(sv);
% for i = 2:length(sv)
%     if sv(i) < tol*sv(1)
%         r = i-1;
%         break
%     end
% end
% 
% U = U(:,1:r);
% [~,R,nflops,r,p] = mgsclpv0(U','tol',0); % use 0 to obtain only ordering
% J = p(1:r); I = p(r+1:end);
% 
% if nargout > 1 % ~istilde(2)
%     E = R(:,J)\R(:,I); % in mgsclpv, cols of R are not actully permuted
%     E = E';
% end

% [J,E,p,nflops] = rrqr2id(U,tol);

% % check error
% U = [eye(r); E]; U = U(invp(p),:);
% V = A(J,:);
% maxerr(U*V,A)/max(max(abs(A)))
