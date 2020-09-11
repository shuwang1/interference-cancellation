function x = TLS(A, y)

[row, col] = size( A );
[U, D, V] = svd( [A y] );
tol = max( size( [A y] ) )*D(1,1)*eps;
K = sum( diag(D) > tol );

x = ( A'*A - D(K,K)^2*eye(col) )\(A'*y);
return

g = V( 1, 1:(K-1) )';
V1 = V( 2:(col+1), 1:(K-1) );
c = V( 1, K:(col+1) )';
V2 = V( 2:(col+1), K:(col+1) );
% x = V1*g/(1-g'*g);
x = V2*c/(c'*c);
return