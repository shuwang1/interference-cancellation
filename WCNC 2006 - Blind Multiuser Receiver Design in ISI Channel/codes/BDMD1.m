%-------------------------------
%
%
%-------------------------------

clear all;
K = 3;

L = 6;

N = 3

a = [2; 8; 8];
A = diag(a);

s1 = [+1; -1; +1; -1; +1; -1];
s1 = s1/norm(s1);
s2 = [+1; +1; +1; -1; +1; -1];
s2 = s2/norm(s2);
s3 = [+1; +1; -1; +1; +1; -1];
s3 = s3/norm(s3);

S = [s1 s2 s3]
S'*S

B = (rand(K, N+1)>0.5)*2 - 1;

sigma = 0;

R0 = zeros(L,N);
R1 = R0;

for n = 1:N
    r = S*A*B(:,n) + sigma*rand(L,1);
    R0(1:L, n) = r;
    R1(1:L, n) = r;
end


r = S*A*B(:,N+1) + sigma*rand(L,1);
R1(1:L, N) = r;

R1\R0