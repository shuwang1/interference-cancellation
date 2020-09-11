%-------------------------------
%
%
%-------------------------------

clear all;
close all;
K = 2;

L = 8;

s1 = [+1; -1; +1; -1; +1; -1; +1; -1];
s1 = s1/norm(s1);
s2 = [+1; +1; +1; -1; +1; -1; +1; -1];
s2 = s2/norm(s2);

S = [s1 s2];
Rs=S'*S
Rs_=inv(S'*S);

sigma=1;

%SNR = 2;
SNR = [0.5; 1; 1.5; 2; 2.5; 3; 3.5; 4];

% A21 = [0.2; 0.4; 0.6; 0.8; 1; 1.5; 2; 2.5; 3; 3.5; 4];
A21 = [0.1; 1; 2];

P = eye(L)-s1*s1';

TotalTrails = 1000;

for n = 1:length(A21)
    for m = 1:length(SNR)
    
        BER_BDD(m,n) = 0;
        BER_DD(m,n) = 0;
        BER_MF(m,n) = 0;
    
        a = [1; A21(n,1)*SNR(m,1)];
        A = diag(a);
        
        w_DD  = (Rs_(1,:)*S')';
        w_MF  = s1;

        for trial = 1: TotalTrails
            
            B = [[1;0] [1;-1]*( (rand > 0.5)*2 -1 )];
            
            S_B  = S*A*B;
            % S_B  = S*A*B + sigma*randn(L,K);

            b = (rand(K,2)>0.5)*2 -1;
            r =S*A*b+sigma*randn(L,2);
            
            d = (inv(S_B'*S_B)*S_B') *r;
            
            dd = d(1,2)/d(1,1);

            x=(1 + B(1,2)*d(2,1) )*dd - B(1,2)*d(2,2);
            x_=(-1 + B(1,2)*d(2,1) )*dd - B(1,2)*d(2,2);
            
            if ( abs(abs(x_)-1) > abs(abs(x)-1) )
                b_BDD = sign(x);
            else
                b_BDD = sign(-x_);
            end
            
            b_DD = sign(w_DD' * r(:,2))/sign(w_DD' * r(:,1));
            b_MF = sign(w_MF' * r(:,2))/sign(w_MF' * r(:,1));
   
            if ( (b_BDD*( b(1,1)/b(1,2) )) < 0 )
                BER_BDD(m,n) = BER_BDD(m,n) + 1;
            end
 
            if ( (b_DD*( b(1,1)/b(1,2) )) < 0 )
                BER_DD(m,n) = BER_DD(m,n) + 1;
            end

            if ( (b_MF*( b(1,1)/b(1,2) )) < 0 )
                BER_MF(m,n) = BER_MF(m,n) + 1;
            end
            

        end
    end
end

BER_MF = BER_MF/TotalTrails
BER_DD = BER_DD/TotalTrails
BER_BDD = BER_LS/TotalTrails
