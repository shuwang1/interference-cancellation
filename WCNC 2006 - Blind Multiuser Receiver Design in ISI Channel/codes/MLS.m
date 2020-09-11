function x = MLS(A, y, G)

[row, col] = size(A);

[Q, R] = qr( A(:,1:G) );

A_qr = Q'*[A y];

x(G+1:col,1) = TLS( A_qr(G+1:row,G+1:col), A_qr(G+1:row, col+1) );
x(1:G,1) = A_qr(1:G,1:G)\(A_qr(1:G, col+1) - A_qr(1:G,G+1:col)*x(G+1:col,1));