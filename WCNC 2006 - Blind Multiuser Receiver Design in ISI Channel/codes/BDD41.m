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

w_DD  = (Rs_(1,:)*S')';
w_MF  = s1;

sigma=1;

%SNR = 2;
SNR = [0.5; 1; 1.5; 2; 2.5; 3; 3.5; 4];

% A21 = [0.2; 0.4; 0.6; 0.8; 1; 1.5; 2; 2.5; 3; 3.5; 4];
A21 = [0.1; 1; 10];

P = eye(L)-s1*s1';

TotalTrails = 1000;

for n = 1:length(A21)
    for m = 1:length(SNR)
    
        BER_LS(m,n) = 0;
        BER_TLS(m,n) = 0;
        BER_DD(m,n) = 0;
        BER_MF(m,n) = 0;
    
        a = [ 1; A21(n,1)];
        A = diag(a);
        
        sigma = 1/SNR(m,1);
        P_DD_Mean(m,n)=0;
        P_MF_Mean(m,n)=0;
        P_LS_Mean(m,n)=0;
        P_TLS_Mean(m,n)=0;
        
        P_DD_Std(m,n)=0;
        P_MF_Std(m,n)=0;
        P_LS_Std(m,n)=0;
        P_TLS_Std(m,n)=0;

        
        for trial = 1: TotalTrails
            
            
            B1 = [[1/A(1,1);0] [1;-1]];
            B2 = [[1/A(1,1);0] [1; 1]];
            
            %S_B1 = S*A*B1;
            %S_B2 = S*A*B2;
            S_B1 = S*A*B1 + [zeros(L,1) sigma*randn(L,K-1)];
            S_B2 = S*A*B2 + [zeros(L,1) sigma*randn(L,K-1)];
            
            b = (rand(K,1)>0.5)*2 -1;
            r =S*A*b + sigma*randn(L,1);
            
            [Us, Ss, Vs]=svd([S_B1 r]);
            s_svd1=sort(diag(Ss));

            [Us, Ss, Vs]=svd([S_B2 r]);
            s_svd2=sort(diag(Ss));
            
            d1_LS = ( inv(S_B1'*S_B1)*S_B1' )*r;
            d2_LS = ( inv(S_B2'*S_B2)*S_B2' )*r;

            d1_TLS = ( pinv( S_B1'*S_B1- s_svd1(1,1)^2*eye(K) )*S_B1' )*r;
            d2_TLS = ( pinv( S_B2'*S_B2- s_svd2(1,1)^2*eye(K) )*S_B2' )*r;
            
            delta_d_LS = d2_LS - d1_LS;
            delta_d_TLS = d2_TLS - d1_TLS;
            
            A1_LS = abs( delta_d_LS(1,1)/(d1_LS(2,1)*B1(1,2)-d2_LS(2,1)*B2(1,2) ) );
            P1_LS(trial,1) =A1_LS^2;
            A1_TLS = abs( delta_d_TLS(1,1)/(d1_TLS(2,1)*B1(1,2)-d2_TLS(2,1)*B2(1,2) ) );
            P1_TLS(trial,1) =A1_TLS^2;
            
            b_LS = sign( [1/A1_LS B1(1,2)]*d1_LS + [1/A1_LS B2(1,2)]*d2_LS );
            b_TLS = sign( [1/A1_TLS B1(1,2)]*d1_TLS + [1/A1_TLS B2(1,2)]*d2_TLS );
            
            
            b_DD = w_DD' * r(:,1);
            P_DD(trial,1) = b_DD^2;
            b_MF = w_MF' * r(:,1);
            P_MF(trial,1) = b_MF^2;
   
            if ( ( b_LS*b(1,1) ) < 0 )
                BER_LS(m,n) = BER_LS(m,n) + 1;
            end

            if ( ( b_TLS*b(1,1) ) < 0 )
                BER_TLS(m,n) = BER_TLS(m,n) + 1;
            end
 
            if ( ( b_DD*b(1,1) ) < 0 )
                BER_DD(m,n) = BER_DD(m,n) + 1;
            end

            if ( ( b_MF*b(1,1) ) < 0 )
                BER_MF(m,n) = BER_MF(m,n) + 1;
            end
            

        end
        
        P_LS_Mean(m,n) = mean(P1_LS);
        P_TLS_Mean(m,n) = mean(P1_TLS);
        P_DD_Mean(m,n) = mean(P_DD);
        P_MF_Mean(m,n) = mean(P_MF);        
        
        P_LS_Std(m,n) = std(P1_LS);
        P_TLS_Std(m,n) = std(P1_TLS);
        P_DD_Std(m,n) = std(P_DD);
        P_MF_Std(m,n) = std(P_MF);        

    end
end

BER_MF = BER_MF/TotalTrails
BER_DD = BER_DD/TotalTrails
BER_LS = BER_LS/TotalTrails
BER_TLS = BER_TLS/TotalTrails

%%---------------------------------------------------------------------------%%
SNR = log10(SNR).*20;

A2A1 = log10(A21).*20;

%%-----------------------------------------------------------------------------
figure(1)
semilogy(SNR, P_MF_Std(:,1),'-*', SNR, P_DD_Std(:,1),'-o',SNR, P_LS_Std(:,1),'-^', SNR, P_TLS_Std(:,1),'-v')
grid;
title('The Perofmance of The Power Estimation Schemes with A_2/A_1=0.1');
xlabel('SNR (dB)');
ylabel('STD of the estimation of P_1');
legend('MF','DD','LS','TLS');


figure(2)
semilogy(SNR, P_MF_Std(:,3),'-*', SNR, P_DD_Std(:,3),'-o',SNR, P_LS_Std(:,3),'-^', SNR, P_TLS_Std(:,3),'-v')
grid;
title('The Perofmance of The Power Estimation Schemes with A_2/A_1=10');
xlabel('SNR (dB)');
ylabel('STD of the estimation of P_1');
legend('MF','DD','LS','TLS');

figure(3)
semilogy(SNR, P_MF_Mean(:,1),'-*', SNR, P_DD_Mean(:,1),'-o',SNR, P_LS_Mean(:,1),'-^', SNR, P_TLS_Mean(:,1),'-v')
grid;
title('The Perofmance of The Power Estimation Schemes with A_2/A_1=0.1');
xlabel('SNR (dB)');
ylabel('Mean of the estimation of P_1');
legend('MF','DD','LS','TLS');



figure(4)
semilogy(SNR, P_MF_Mean(:,3),'-*', SNR, P_DD_Mean(:,3),'-o',SNR, P_LS_Mean(:,3),'-^', SNR, P_TLS_Mean(:,3),'-v')
grid;
title('The Perofmance of The Power Estimation Schemes with A_2/A_1=10');
xlabel('SNR (dB)');
ylabel('Mean of the estimation of P_1');
legend('MF','DD','LS','TLS');


%%------------------------------------------------------------------------------------
figure(5)
semilogy(SNR, BER_MF(:,1),'-*',SNR, BER_MF(:,2),'-+',SNR, BER_MF(:,3),'-x', SNR, BER_DD(:,3),'-o', SNR, BER_LS(:,1),'-d', SNR, BER_LS(:,2),'-p', SNR, BER_LS(:,3),'-h')
ylim([1.0/TotalTrails 0.55])
xlim([min(SNR) max(SNR)])
grid;
title('BER of Different Multiuser Detection Algorithms');
xlabel('SNR (dB)');
ylabel('Bit-Error-Rate BER');
legend('MF with A_2/A_1=0.1','MF with A_2/A_1=1','MF with A_2/A_1=2','Decorr. Detector','LS with A_2/A_1=0.1','LS with A_2/A_1=1','LS with A_2/A_1=2');

figure(6)
semilogy(SNR, BER_MF(:,1),'-*',SNR, BER_MF(:,2),'-+',SNR, BER_MF(:,3),'-x', SNR, BER_DD(:,3),'-o', SNR, BER_TLS(:,1),'-d', SNR, BER_TLS(:,2),'-p', SNR, BER_TLS(:,3),'-h')
ylim([1.0/TotalTrails 0.55])
xlim([min(SNR) max(SNR)])
grid;
title('BER of Different Multiuser Detection Algorithms');
xlabel('SNR (dB)');
ylabel('Bit-Error-Rate BER');
legend('MF with A_2/A_1=0.1','MF with A_2/A_1=1','MF with A_2/A_1=2','Decorr. Detector','TLS with A_2/A_1=0.1','TLS with A_2/A_1=1','TLS with A_2/A_1=2');