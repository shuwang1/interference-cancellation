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

TotalTrails = 10000;

for n = 1:length(A21)
    for m = 1:length(SNR)
    
        BER_LS(m,n)  = 0;
        BER_TLS(m,n) = 0;
        BER_DD(m,n) = 0;
        BER_MF(m,n) = 0;
        BER_LS2(m,n)=0; % LS-MMSE 
        BER_test(m,n)=0; % TLS-MMSE 
    
        a = [1; A21(n,1)]*SNR(m,1);
        A = diag(a);
        
        w_DD  = Rs_(1,:)*S';
        w_MF  = s1';

        for trial = 1: TotalTrails
            
            BB = [[1;1]  [1;-1] ]*( (rand(1)>0.5)*2 -1 );
            
            % S_BB  = S*A*BB;
            S_BB  = S*A*BB + sigma*randn(L,K);
            S_BB_ = P*S_BB;

            b = (rand(K,1)>0.5)*2 -1;
            r =S*A*b+sigma*randn(L,1);
            r_ = P*r;
                        
            
            [Us, Ss, Vs]=svd([S_BB_ r_]);
            s_svd=sort(diag(Ss));
            
            [Us_, Ss_, Vs_]=svd(S_BB_);
            
            error=s_svd(1,1)^2;
            
            for p =1:K
                error = error + s_svd(1,1)^2*(Us(:,p)'*r)^2/(Ss_(p,p)^2-s_svd(1,1)^2);
            end
            
        
            b_LS  = sign( s1'*r - s1'*S_BB*( pinv(S_BB_)*r_ ) );
            b_TLS = sign( s1'*r - s1'*S_BB*pinv( S_BB_'*S_BB_ - s_svd(1,1)^2*eye(K) )*S_BB_'*r_ );
            b_DD = sign(w_DD*r);
            b_MF = sign(w_MF*r);
   
            if ( (b_LS*b(1,1)) < 0 )
                BER_LS(m,n) = BER_LS(m,n) + 1;
            end
 
            if ( (b_TLS*b(1,1)) < 0 )
                BER_TLS(m,n) = BER_TLS(m,n) + 1;
            end

            if ( (b_DD*b(1,1)) < 0 )
                BER_DD(m,n) = BER_DD(m,n) + 1;
            end

            if ( (b_MF*b(1,1)) < 0 )
                BER_MF(m,n) = BER_MF(m,n) + 1;
            end
            

        end
    end
end

BER_MF = BER_MF/TotalTrails;
BER_DD = BER_DD/TotalTrails;
BER_LS = BER_LS/TotalTrails;
BER_TLS = BER_TLS/TotalTrails;

%%------------------------------------------------------------------------------------------------%%

figure(1)
loglog(SNR, BER_MF(:,1),'-*', SNR, BER_DD(:,1),'-o', SNR, BER_LS(:,1),'-+', SNR, BER_TLS(:,1),'-d')
ylim([1.0/TotalTrails 0.5])
xlim([min(SNR) max(SNR)])
grid;
title('BER of Different Multiuser Detection Algorithms');
xlabel('A_k/\sigma_n');
ylabel('Bit-Error-Rate BER');
legend('MF','DD','LS','TLS');

figure(2)
loglog(SNR, BER_MF(:,2),'-*', SNR, BER_DD(:,2),'-o', SNR, BER_LS(:,2),'-+', SNR, BER_TLS(:,2),'-d')
ylim([1.0/TotalTrails 0.5])
xlim([min(SNR) max(SNR)])
grid;
title('BER of Different Multiuser Detection Algorithms');
xlabel('A_k/\sigma_n');
ylabel('Bit-Error-Rate BER');
legend('MF','DD','LS','TLS');

figure(3)
loglog(SNR, BER_MF(:,3),'-*', SNR, BER_DD(:,3),'-o', SNR, BER_LS(:,3),'-+', SNR, BER_TLS(:,3),'-d')
ylim([1.0/TotalTrails 0.5])
xlim([min(SNR) max(SNR)])
grid;
title('BER of Different Multiuser Detection Algorithms');
xlabel('A_k/\sigma_n');
ylabel('Bit-Error-Rate BER');
legend('MF','DD','LS','TLS');

SNR = log10(SNR).*20;

figure(4)
semilogy(SNR, BER_MF(:,1),'-*',SNR, BER_MF(:,2),'-+',SNR, BER_MF(:,3),'-x', SNR, BER_DD(:,3),'-o', SNR, BER_LS(:,3),'-s', SNR, BER_TLS(:,1),'-d', SNR, BER_TLS(:,2),'-p', SNR, BER_TLS(:,3),'-h')
ylim([1.0/TotalTrails 0.55])
xlim([min(SNR) max(SNR)])
grid;
title('BER of Different Multiuser Detection Algorithms');
xlabel('SNR (dB)');
ylabel('Bit-Error-Rate BER');
legend('MF with A_2/A_1=0.1','MF with A_2/A_1=1','MF with A_2/A_1=2','Decorrelator Detector','LS Semi-Blind DD','TLS with A_2/A_1=0.1','TLS with A_2/A_1=1','TLS with A_2/A_1=2');

figure(5)
semilogy(SNR, BER_MF(:,1),'-*',SNR, BER_MF(:,2),'-+',SNR, BER_MF(:,3),'-x', SNR, BER_DD(:,3),'-o', SNR, BER_LS(:,1),'-d', SNR, BER_LS(:,2),'-p', SNR, BER_LS(:,3),'-h')
ylim([1.0/TotalTrails 0.55])
xlim([min(SNR) max(SNR)])
grid;
title('BER of Different Multiuser Detection Algorithms');
xlabel('SNR (dB)');
ylabel('Bit-Error-Rate BER');
legend('MF with A_2/A_1=0.1','MF with A_2/A_1=1','MF with A_2/A_1=2','Decorrelator Detector','LS with A_2/A_1=0.1','LS with A_2/A_1=1','LS with A_2/A_1=2');

figure(6)
semilogy(SNR, BER_MF(:,1),'-*',SNR, BER_MF(:,2),'-+',SNR, BER_MF(:,3),'-x', SNR, BER_DD(:,3),'-o', SNR, BER_TLS(:,1),'-d', SNR, BER_TLS(:,2),'-p', SNR, BER_TLS(:,3),'-h')
ylim([1.0/TotalTrails 0.55])
xlim([min(SNR) max(SNR)])
grid;
title('BER of Different Multiuser Detection Algorithms');
xlabel('SNR (dB)');
ylabel('Bit-Error-Rate BER');
legend('MF with A_2/A_1=0.1','MF with A_2/A_1=1','MF with A_2/A_1=2','Decorrelator Detector','TLS with A_2/A_1=0.1','TLS with A_2/A_1=1','TLS with A_2/A_1=2');