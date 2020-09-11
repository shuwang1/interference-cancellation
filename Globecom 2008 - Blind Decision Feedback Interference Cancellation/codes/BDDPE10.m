%-------------------------------------------------------------------
%   
%   One-Shot Blind Decorrelating Detection and Power Estimation
% for Synchronous CDMA
%
%-------------------------------------------------------------------
%   BDDPE10.m
%
%   Created on 05/25/2002 
%       ---- the d estimation, power estimation, BER and NFR 
%
%-------------------------------------------------------------------

clear all;
close all;

K = 2;
M=2;

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

SNR = 10^(10/20);
A21 = [0.1; 0.2; 0.4; 0.6; 0.8; 1; 2; 4; 6; 8; 10];

%SNR = [0.1; 0.2; 0.4; 0.6; 0.8; 1; 2; 3; 4; 5];
%A21 = [0.1; 1; 10];

P = eye(L)-s1*s1'/norm(s1);
[Q, R] = qr(s1)
sigma=1;

TotalTrails = 20000

for n = 1:length(A21)
    for m = 1:length(SNR)

        BER_DD(m,n) = 0;
        BER_MF(m,n) = 0;
        BER_LS(m,n) = 0;
        BER_TLS(m,n) = 0;
        BER_MLS(m,n) = 0;        

        P_DD_Mean(m,n)=0;
        P_MF_Mean(m,n)=0;
        P_LS_Mean(m,n)=0;
        P_TLS_Mean(m,n)=0;
        P_MLS_Mean(m,n)=0;
        
        P_DD_Std(m,n)=0;
        P_MF_Std(m,n)=0;
        P_LS_Std(m,n)=0;
        P_TLS_Std(m,n)=0;
        P_MLS_Std(m,n)=0;
        
        a = [ 1; A21(n,1)];
        A = diag(a);
        
        sigma = 1/SNR(m,1);

        trial = 1;
        
        while trial <= TotalTrails

            b = (rand(K,1)>0.5)*2 -1;
            r = S*A*b + sigma*randn(L,1);
            
            B1 = [ [1/A(1,1);0] [1; 1]*( (rand > 0.5)*2 -1 ) ];
            B2 = [ [1/A(1,1);0] [1;-1]*( (rand > 0.5)*2 -1 ) ];
            
            d1 = pinv(B1)*b;
            d2 = pinv(B2)*b;

            if ( abs( d1(2,1)*B1(1,2)-d2(2,1)*B2(1,2) ) < 0.01 )
                continue
            end

            %S_B1 = S*A*B1;
            %S_B2 = S*A*B2;
            
            S_B1 = S*A*B1 + [zeros(L,1) sigma*randn(L,M-1)];
            S_B2 = S*A*B2 + [zeros(L,1) sigma*randn(L,M-1)];
            
            %%===============================================%%
            
            d1_LS = ( inv(S_B1'*S_B1)*S_B1' )*r;
            d2_LS = ( inv(S_B2'*S_B2)*S_B2' )*r;

            dd_LS(2*trial-1,1)=norm(d1_LS-d1);
            dd_LS(2*trial,1)=norm(d2_LS-d2);
            
            delta_d_LS  = d2_LS  - d1_LS;
            A1_LS(trial,1)  = abs( delta_d_LS(1,1)/( B1(1,2:M)*d1_LS(2:M,1)-B2(1,2:M)*d2_LS(2:M,1) ) );
            P1_LS(trial,1)  = A1_LS(trial,1)^2;
            if ( P1_LS(trial,1)>1000 )
                continue
            end

            b_LS = sign( [1/A1_LS(trial,1) B1(1,2:M)]*d1_LS + [1/A1_LS(trial,1) B2(1,2:M)]*d2_LS );
            
            %================================================%
            [Us1, Ss1, Vs1]=svd([S_B1 r]);
            s_svd1=sort(diag(Ss1));
            d1_TLS = ( pinv( S_B1'*S_B1 - (s_svd1(1,1)^2)*eye(M) )*S_B1' )*r;
            
            [Us2, Ss2, Vs2]=svd([S_B2 r]);
            s_svd2=sort(diag(Ss2));
            d2_TLS = ( pinv( S_B2'*S_B2 - (s_svd2(1,1)^2)*eye(M) )*S_B2' )*r;
            
            dd_TLS(2*trial-1,1) = norm(d1_TLS-d1);            
            dd_TLS(2*trial  ,1) = norm(d2_TLS-d2);            
            
            delta_d_TLS = d2_TLS - d1_TLS;
            A1_TLS(trial,1) = abs( delta_d_TLS(1,1)/( B1(1,2:M)*d1_TLS(2:M,1)-B2(1,2:M)*d2_TLS(2:M,1) ) );
            P1_TLS(trial,1) = A1_TLS(trial,1)^2;
            
            if ( P1_TLS(trial,1)>1000 )
                continue
            end
                
            b_TLS = sign( [1/A1_TLS(trial,1) B1(1,2:M)]*d1_TLS + [1/A1_TLS(trial,1) B2(1,2:M)]*d2_TLS );
            

            %=================================================%
            
            temp = Q'*[S_B1 r];
            [Us1, Ss1, Vs1]=svd(temp(2:L,2:M+1));
            s_svd1=sort(diag(Ss1));
            d1_MLS = ( pinv( S_B1'*S_B1 - (s_svd1(1,1)^2)*diag([0 ones(1,M-1)]) )*S_B1' )*r;
            
            temp = Q'*[S_B2 r];
            [Us2, Ss2, Vs2]=svd(temp(2:L,2:M+1));
            s_svd2=sort(diag(Ss2));
            d2_MLS = ( pinv( S_B2'*S_B2 - (s_svd2(1,1)^2)*diag([0 ones(1,M-1)]) )*S_B2' )*r;

            dd_MLS(2*trial-1,1) = norm(d1_MLS-d1);            
            dd_MLS(2*trial  ,1) = norm(d2_MLS-d2);            
            
            delta_d_MLS = d2_MLS - d1_MLS;
            A1_MLS(trial,1) = abs( delta_d_MLS(1,1)/( B1(1,2:M)*d1_MLS(2:M,1)-B2(1,2:M)*d2_MLS(2:M,1) ) );
            P1_MLS(trial,1) = A1_MLS(trial,1)^2;
            
            if ( P1_MLS(trial,1)>1000 )
                continue
            end

            b_MLS = sign( [1/A1_MLS(trial,1) B1(1,2:M)]*d1_MLS + [1/A1_MLS(trial,1) B2(1,2:M)]*d2_MLS );

            %=================================================%
                
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

            if ( ( b_MLS*b(1,1) ) < 0 )
                BER_MLS(m,n) = BER_MLS(m,n) + 1;
            end

            
            if ( ( b_DD*b(1,1) ) < 0 )
                BER_DD(m,n) = BER_DD(m,n) + 1;
            end
                
            if ( ( b_MF*b(1,1) ) < 0 )
                BER_MF(m,n) = BER_MF(m,n) + 1;
            end
            
            trial = trial + 1;

        end
        
        MEAN_dd_LS(n,m)=mean(dd_LS);
        MEAN_dd_TLS(n,m)=mean(dd_TLS);
        MEAN_dd_MLS(n,m)=mean(dd_MLS);

        STD_dd_LS(n,m)=std(dd_LS);
        STD_dd_TLS(n,m)=std(dd_TLS);
        STD_dd_MLS(n,m)=std(dd_MLS);

        P_LS_Mean(m,n) = mean(P1_LS);
        P_TLS_Mean(m,n) = mean(P1_TLS);
        P_MLS_Mean(m,n) = mean(P1_MLS);
        P_DD_Mean(m,n) = mean(P_DD);
        P_MF_Mean(m,n) = mean(P_MF);        
        
        P_LS_Std(m,n) = std(P1_LS);
        P_TLS_Std(m,n) = std(P1_TLS);
        P_MLS_Std(m,n) = std(P1_MLS);        
        P_DD_Std(m,n) = std(P_DD);
        P_MF_Std(m,n) = std(P_MF);        

    end
end

BER_MF = BER_MF/TotalTrails
BER_DD = BER_DD/TotalTrails
BER_LS = BER_LS/TotalTrails
BER_TLS = BER_TLS/TotalTrails
BER_MLS = BER_MLS/TotalTrails

%%---------------------------------------------------------------------------%%
P1_SNR = SNR.^2;
SNR = log10(SNR).*20;

A2A1 = log10(A21).*20;

%%---------------------------------------------------------------
figure(1)
semilogy(SNR, STD_dd_LS(1,:),'-<', SNR, STD_dd_TLS(1,:),'->',SNR, STD_dd_LS(3,:),'-v', SNR, STD_dd_TLS(3,:),'-^')
grid;
title('The Perofmance of The Estimation Schemes for d_n');
xlabel('SNR (dB)');
ylabel('STD of  the difference between the estimated d_n and the real d_n');
legend('LS with A_2/A_1=0.1','TLS with A_2/A_1=0.1','LS with A_2/A_1=10','TLS with A_2/A_1=10');

figure(2)
semilogy(SNR, MEAN_dd_LS(1,:),'-<', SNR, MEAN_dd_TLS(1,:),'->',SNR, MEAN_dd_LS(3,:),'-v', SNR, MEAN_dd_TLS(3,:),'-^')
grid;
title('The Perofmance of The Estimation Schemes for d_n');
xlabel('SNR (dB)');
ylabel('Mean of  the difference between the estimated d_n and the real d_n');
legend('LS with A_2/A_1=0.1','TLS with A_2/A_1=0.1','LS with A_2/A_1=10','TLS with A_2/A_1=10');

figure(3)
semilogy(SNR, STD_dd_LS(1,:),'-<', SNR, STD_dd_MLS(1,:),'->',SNR, STD_dd_LS(3,:),'-v', SNR, STD_dd_MLS(3,:),'-^')
grid;
title('The Perofmance of The Estimation Schemes for d_n');
xlabel('SNR (dB)');
ylabel('STD of  the difference between the estimated d_n and the real d_n');
legend('LS with A_2/A_1=0.1','MLS with A_2/A_1=0.1','LS with A_2/A_1=10','TLS with A_2/A_1=10');

figure(4)
semilogy(SNR, MEAN_dd_LS(1,:),'-<', SNR, MEAN_dd_MLS(1,:),'->',SNR, MEAN_dd_LS(3,:),'-v', SNR, MEAN_dd_MLS(3,:),'-^')
grid;
title('The Perofmance of The Estimation Schemes for d_n');
xlabel('SNR (dB)');
ylabel('Mean of  the difference between the estimated d_n and the real d_n');
legend('LS with A_2/A_1=0.1','MLS with A_2/A_1=0.1','LS with A_2/A_1=10','TLS with A_2/A_1=10');



%%-----------------------------------------------------------------------------

figure(5)
semilogy(SNR, P_MF_Std(:,1),'-s', SNR, P_DD_Std(:,1),'-o',SNR, P_LS_Std(:,1),'-^', SNR, P_TLS_Std(:,1),'-v')
grid;
title('The Perofmance of The Power Estimation Schemes with A_2/A_1=0.1');
xlabel('SNR (dB)');
ylabel('STD of the estimation of P_1');
legend('MF','DD','LS','TLS');

figure(6)
semilogy(SNR, P_MF_Std(:,3),'-s', SNR, P_DD_Std(:,3),'-o',SNR, P_LS_Std(:,3),'-^', SNR, P_TLS_Std(:,3),'-v')
grid;
title('The Perofmance of The Power Estimation Schemes with A_2/A_1=10');
xlabel('SNR (dB)');
ylabel('STD of the estimation of P_1');
legend('MF','DD','LS','TLS');

figure(7)
semilogy(SNR, P_MF_Mean(:,1),'-s', SNR, P_DD_Mean(:,1),'-o',SNR, P_LS_Mean(:,1),'-^', SNR, P_TLS_Mean(:,1),'-v')
grid;
title('The Perofmance of The Power Estimation Schemes with A_2/A_1=0.1');
xlabel('SNR (dB)');
ylabel('Mean of the estimation of P_1');
legend('MF','DD','LS','TLS');

figure(8)
semilogy(SNR, P_MF_Mean(:,3),'-s', SNR, P_DD_Mean(:,3),'-o',SNR, P_LS_Mean(:,3),'-^', SNR, P_TLS_Mean(:,3),'-v')
grid;
title('The Perofmance of The Power Estimation Schemes with A_2/A_1=10');
xlabel('SNR (dB)');
ylabel('Mean of the estimation of P_1');
legend('MF','DD','LS','TLS');

figure(9)
semilogy(SNR, P_MF_Std(:,1),'-s', SNR, P_DD_Std(:,1),'-o',SNR, P_LS_Std(:,1),'-^', SNR, P_MLS_Std(:,1),'-v')
grid;
title('The Perofmance of The Power Estimation Schemes with A_2/A_1=0.1');
xlabel('SNR (dB)');
ylabel('STD of the estimation of P_1');
legend('MF','DD','LS','MLS');

figure(10)
semilogy(SNR, P_MF_Std(:,3),'-s', SNR, P_DD_Std(:,3),'-o',SNR, P_LS_Std(:,3),'-^', SNR, P_MLS_Std(:,3),'-v')
grid;
title('The Perofmance of The Power Estimation Schemes with A_2/A_1=10');
xlabel('SNR (dB)');
ylabel('STD of the estimation of P_1');
legend('MF','DD','LS','MLS');

figure(11)
semilogy(SNR, P_MF_Mean(:,1),'-s', SNR, P_DD_Mean(:,1),'-o',SNR, P_LS_Mean(:,1),'-^', SNR, P_MLS_Mean(:,1),'-v')
grid;
title('The Perofmance of The Power Estimation Schemes with A_2/A_1=0.1');
xlabel('SNR (dB)');
ylabel('Mean of the estimation of P_1');
legend('MF','DD','LS','MLS');

figure(12)
semilogy(SNR, P_MF_Mean(:,3),'-s', SNR, P_DD_Mean(:,3),'-o',SNR, P_LS_Mean(:,3),'-^', SNR, P_MLS_Mean(:,3),'-v')
grid;
title('The Perofmance of The Power Estimation Schemes with A_2/A_1=10');
xlabel('SNR (dB)');
ylabel('Mean of the estimation of P_1');
legend('MF','DD','LS','MLS');


figure(13)
semilogy(SNR, P_MF_Std(:,1),'-+', SNR, P_MF_Std(:,3),'-*', SNR, P_DD_Std(:,1),'-o', SNR, P_DD_Std(:,3),'-s',SNR, P_LS_Std(:,1),'-^',SNR, P_LS_Std(:,3),'-v',SNR, P_MLS_Std(:,1),'-<',SNR, P_MLS_Std(:,3),'->')
grid;
title('The Perofmance of The Power Estimation Schemes with A_2/A_1=0.1');
xlabel('SNR (dB)');
ylabel('STD of the estimation of P_1');
legend('MF with A_2/A_1=0.1','MF with A_2/A_1=10','DD with A_2/A_1=0.1','DD with A_2/A_1=10','LS with A_2/A_1=0.1','LS with A_2/A_1=10','MLS with A_2/A_1=0.1','MLS with A_2/A_1=10');

figure(14)
semilogy(SNR,ones(size(SNR)),'-.',SNR, P_MF_Mean(:,1),'-+',SNR, P_MF_Mean(:,3),'-*', SNR, P_DD_Mean(:,1),'-o', SNR, P_DD_Mean(:,3),'-s',SNR, P_LS_Mean(:,1),'-^',SNR, P_LS_Mean(:,3),'-v',SNR, P_MLS_Mean(:,1),'-<',SNR, P_MLS_Mean(:,3),'->')
grid;
title('The Perofmance of The Power Estimation Schemes with A_2/A_1=0.1');
xlabel('SNR (dB)');
ylabel('Mean of the estimation of P_1');
legend('The orginal power','MF with A_2/A_1=0.1','MF with A_2/A_1=10','DD with A_2/A_1=0.1','DD with A_2/A_1=10','LS with A_2/A_1=0.1','LS with A_2/A_1=10','MLS with A_2/A_1=0.1','MLS with A_2/A_1=10');


figure(15)
semilogy(SNR, P_MF_Std(:,1),'-+', SNR, P_MF_Std(:,2),'-*', SNR, P_MF_Std(:,3),'-x', SNR, P_DD_Std(:,1),'-d', SNR, P_DD_Std(:,2),'-o', SNR, P_DD_Std(:,3),'-s',SNR, P_LS_Std(:,1),'-^',SNR, P_LS_Std(:,2),'-v',SNR, P_LS_Std(:,3),'-<')
grid;
title('The Perofmance of The Power Estimation Schemes with A_2/A_1=0.1');
xlabel('SNR (dB)');
ylabel('STD of the estimation of P_1');
legend('MF with A_2/A_1=0.1','MF with A_2/A_1=1.0','MF with A_2/A_1=10','DD with A_2/A_1=0.1','DD with A_2/A_1=1.0','DD with A_2/A_1=10','LS with A_2/A_1=0.1','LS with A_2/A_1=1.0','LS with A_2/A_1=10');

figure(16)
semilogy(SNR,ones(size(SNR)),'-.',SNR, P_MF_Mean(:,1),'-+',SNR, P_MF_Mean(:,2),'-x',SNR, P_MF_Mean(:,3),'-*',SNR, P_DD_Mean(:,1),'-d', SNR, P_DD_Mean(:,2),'-o', SNR, P_DD_Mean(:,3),'-s',SNR, P_LS_Mean(:,1),'-^',SNR, P_LS_Mean(:,2),'-v',SNR, P_LS_Mean(:,3),'-<')
grid;
title('The Perofmance of The Power Estimation Schemes with A_2/A_1=0.1');
xlabel('SNR (dB)');
ylabel('Mean of the estimation of P_1');
legend('The orginal power','MF with A_2/A_1=0.1','MF with A_2/A_1=1.0','MF with A_2/A_1=10','DD with A_2/A_1=0.1','DD with A_2/A_1=1.0','DD with A_2/A_1=10','LS with A_2/A_1=0.1','LS with A_2/A_1=1.0','LS with A_2/A_1=10');


%%------------------------------------------------------------------------------------

figure(17)
semilogy(SNR, BER_MF(:,1),'-*',SNR, BER_MF(:,2),'-+',SNR, BER_MF(:,3),'-x', SNR, BER_DD(:,2),'-o', SNR, BER_LS(:,1),'-d', SNR, BER_LS(:,2),'-p', SNR, BER_LS(:,3),'-h')
ylim([1.0/TotalTrails 0.55])
%xlim([min(SNR) max(SNR)])
grid;
title('BER of Different Multiuser Detection Algorithms');
xlabel('SNR (dB)');
ylabel('Bit-Error-Rate BER');
legend('MF with A_2/A_1=0.1','MF with A_2/A_1=1.0','MF with A_2/A_1=10','DD with A_2/A_1=1.0','LS with A_2/A_1=0.1','LS with A_2/A_1=1.0','LS with A_2/A_1=10');

figure(18)
semilogy(SNR, BER_MF(:,1),'-*',SNR, BER_MF(:,2),'-+',SNR, BER_MF(:,3),'-x', SNR, BER_DD(:,2),'-o', SNR, BER_TLS(:,1),'-d', SNR, BER_TLS(:,2),'-p', SNR, BER_TLS(:,3),'-h')
ylim([1.0/TotalTrails 0.55])
%xlim([min(SNR) max(SNR)])
grid;
title('BER of Different Multiuser Detection Algorithms');
xlabel('SNR (dB)');
ylabel('Bit-Error-Rate BER');
legend('MF with A_2/A_1=0.1','MF with A_2/A_1=1.0','MF with A_2/A_1=10','DD with A_2/A_1=1.0','TLS with A_2/A_1=0.1','TLS with A_2/A_1=1.0','TLS with A_2/A_1=10');

figure(19)
semilogy(SNR, BER_MF(:,1),'-*',SNR, BER_MF(:,2),'-+',SNR, BER_MF(:,3),'-x', SNR, BER_DD(:,2),'-o', SNR, BER_MLS(:,1),'-d', SNR, BER_MLS(:,2),'-p', SNR, BER_MLS(:,3),'-h')
ylim([1.0/TotalTrails 0.55])
%xlim([min(SNR) max(SNR)])
grid;
title('BER of Different Multiuser Detection Algorithms');
xlabel('SNR (dB)');
ylabel('Bit-Error-Rate BER');
legend('MF with A_2/A_1=0.1','MF with A_2/A_1=1.0','MF with A_2/A_1=10','DD with A_2/A_1=1.0','MLS with A_2/A_1=0.1','MLS with A_2/A_1=1.0','MLS with A_2/A_1=10');


%%------------------------------------------------------------------------------------------
figure(20)
semilogy(A2A1, BER_MF(1,:),'-*', A2A1, BER_DD(1,:),'-o', A2A1, BER_LS(1,:),'-^', A2A1, BER_TLS(1,:),'-v')
grid;
title('Near-Far Resistance Performance of Different Detectors');
xlabel('Near-Far Resistance A_2/A_1(dB)');
ylabel('Bit-Error-Rate (BER)');
legend('MF','DD','LS','TLS');

figure(21)
semilogy(A2A1, BER_MF(1,:),'-*', A2A1, BER_DD(1,:),'-o', A2A1, BER_LS(1,:),'-^', A2A1, BER_MLS(1,:),'-v')
grid;
title('Near-Far Resistance Performance of Different Detectors');
xlabel('Near-Far Resistance A_2/A_1(dB)');
ylabel('Bit-Error-Rate (BER)');
legend('MF','DD','LS','MLS');

figure(22)
semilogy(A2A1, BER_MF(1,:),'-*', A2A1, BER_DD(1,:),'-o', A2A1, BER_LS(1,:),'-^')
grid;
title('Near-Far Resistance Performance of Different Detectors');
xlabel('Near-Far Resistance A_2/A_1(dB)');
ylabel('Bit-Error-Rate (BER)');
legend('MF','DD','LS');

%%==================================================================================================================%
figure(23)
semilogy(A2A1,ones(size(A2A1)),'-.',A2A1, P_MF_Mean(1,:),'-*', A2A1, P_DD_Mean(1,:),'-o', A2A1, P_LS_Mean(1,:),'-^', A2A1, P_MLS_Mean(1,:),'-v')
grid;
title('Near-Far Resistance Performance of Different Detectors');
xlabel('Near-Far Resistance A_2/A_1(dB)');
ylabel('Mean of the estimation of P_1');
legend('The original power','MF','DD','LS','MLS');

figure(24)
semilogy(A2A1, P_MF_Std(1,:),'-*', A2A1, P_DD_Std(1,:),'-o', A2A1, P_LS_Std(1,:),'-^', A2A1, P_MLS_Std(1,:),'-v')
grid;
title('Near-Far Resistance Performance of Different Detectors');
xlabel('Near-Far Resistance A_2/A_1(dB)');
ylabel('STD of the estimation of P_1');
legend('MF','DD','LS','MLS');


figure(25)
semilogy(A2A1,ones(size(A2A1)),'-.',A2A1, P_MF_Mean(1,:),'-*', A2A1, P_DD_Mean(1,:),'-o', A2A1, P_LS_Mean(1,:),'-^')
grid;
title('Near-Far Resistance Performance of Different Detectors');
xlabel('Near-Far Resistance A_2/A_1(dB)');
ylabel('Mean of the estimation of P_1');
legend('The original power','MF','DD','LS');

figure(26)
semilogy(A2A1, P_MF_Std(1,:),'-*', A2A1, P_DD_Std(1,:),'-o', A2A1, P_LS_Std(1,:),'-^')
grid;
title('Near-Far Resistance Performance of Different Detectors');
xlabel('Near-Far Resistance A_2/A_1(dB)');
ylabel('STD of the estimation of P_1');
legend('MF','DD','LS');

%%=============================================================================================
close all;
figure(31)
semilogy(A2A1, BER_MF(1,:),'-', A2A1, BER_DD(1,:),'--', A2A1, BER_LS(1,:),'-s', A2A1, BER_TLS(1,:),'-d', A2A1, BER_MLS(1,:),'-o')
grid;
title('Near-Far Resistance Performance of Different Detectors');
xlabel('Near-Far Resistance A_2/A_1(dB)','FontWeight','bold','FontSize',14);
ylabel('BER','FontWeight','bold','FontSize',14);
legend('MF','DD','LS','TLS','MLS');
