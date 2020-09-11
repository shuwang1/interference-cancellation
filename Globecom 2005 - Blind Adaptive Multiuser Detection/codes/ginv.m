function B = ginv(A, K)
%   general inverse
%
%   
%
%   Started by Shu Wang, 12/09/2004



[U, S, V] = svd(A);
B = U(1:,1:K)*S(1:K,1:K)*V(:,1:K)';