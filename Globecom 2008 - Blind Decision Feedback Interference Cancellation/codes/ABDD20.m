
%------------------------------------------------------------------
% One-Shot Blind decorrelating detector for Asynchronous CDMA
% 
% Created from ASBD15.m on 08/17/2003
%
%
%------------------------------------------------------------------

clear all;
close all;

K = 2;
L = 24;
P = 1;
M = (P+1)*K-1;

s1 = [+1; -1; +1; -1; +1; -1; +1; -1; +1; -1; +1; -1; +1; -1; +1; -1; +1; -1; +1; -1; +1; -1; +1; -1];
s1 = s1/sqrt(L);
s2 = [+1; +1; +1; -1; +1; -1; +1; -1; +1; +1; +1; -1; +1; +1; +1; -1; +1; -1; +1; -1; +1; +1; +1; -1];
s2 = s2/sqrt(L);
S = [s1 s2];
S'*S;

s21 = [-1; +1; -1; +1; -1; +1; +1; +1; -1; +1; +1; +1; -1; +1; -1; +1; -1; +1; +1; +1; -1;  0;0;0 ];
s21 = s21/sqrt(L);
s22 = [ zeros(21,1);  +1; +1; +1];
s22 = s22/sqrt(L);

SS_1 = [s1 s21 s22];
SS_1'*SS_1;
R_SS_1_ = inv(SS_1'*SS_1);

SS = zeros(P*L, (P+1)*K-1);
SS(:,1:P) = kron(eye(P),s1);

tao=[0 3];

for k =2:K
    SS(1:L-tao(k),(k-1)*(P+1) ) = S(tao(k)+1:L,k);
    % SS(L-tao(k)+1:L-tao(k)+(P-1)*L, (k-1)*(P+1)+1:(k-1)*(P+1)+P-1) = kron(eye(P-1), S(:,k));
    SS(L-tao(k)+(P-1)*L+1:P*L, (k-1)*(P+1)+P) = S(1:tao(k),k);
end

sigma=1;

SNR = [0.1; 0.2; 0.4; 0.8; 1; 2; 3; 4; 5];
A21 = [0.1; 1; 10];

%SNR = [2];
%A21 = [0.1; 0.2; 0.4; 0.6; 0.8; 1; 1.5; 2; 2.5; 3; 3.5; 4; 5];

[Q,R]=qr(SS(:,1:P));
Q'*SS;

TotalTrails = 200;

w_DD = (R_SS_1_(1,:)*SS_1');
w_MF = s1';

for m = 1:length(SNR)
    for n = 1:length(A21)

        BER_DD1(m,n)  = 0;
        BER_MF1(m,n) = 0;

        BER_DD(m,n)  = 0;
        BER_MF(m,n) = 0;
        BER_LS(m,n) = 0;
        BER_TLS(m,n) = 0;
        BER_MLS(m,n) = 0;        
        
        
        BER_LS_LS1(m,n)  = 0;
        BER_LS_LS2(m,n)  = 0;
        BER_LS_TLS(m,n)  = 0;

        BER_TLS_LS1(m,n)  = 0;
        BER_TLS_LS2(m,n)  = 0;
        BER_TLS_TLS(m,n)  = 0;

        BER_MLS_LS1(m,n)  = 0;
        BER_MLS_LS2(m,n)  = 0;
        BER_MLS_TLS(m,n)  = 0;
        
        P_DD_Mean(m,n)=0;
        P_MF_Mean(m,n)=0;
        
        P_LS_LS1_Mean(m,n)=0;
        P_LS_LS2_Mean(m,n)=0;
        P_LS_TLS_Mean(m,n)=0;
        P_TLS_LS1_Mean(m,n)=0;
        P_TLS_LS2_Mean(m,n)=0;
        P_TLS_TLS_Mean(m,n)=0;
        P_MLS_LS1_Mean(m,n)=0;
        P_MLS_LS2_Mean(m,n)=0;
        P_MLS_TLS_Mean(m,n)=0;
        
        P_DD_Std(m,n)=0;
        P_MF_Std(m,n)=0;
        P_LS_LS1_Std(m,n)=0;
        P_LS_LS2_Std(m,n)=0;
        P_LS_TLS_Std(m,n)=0;
        P_TLS_LS1_Std(m,n)=0;
        P_TLS_LS2_Std(m,n)=0;
        P_TLS_TLS_Std(m,n)=0;
        P_MLS_LS1_Std(m,n)=0;
        P_MLS_LS2_Std(m,n)=0;
        P_MLS_TLS_Std(m,n)=0;
    
        a = [ones(P,1); A21(n,1)*ones(P+1,1)]*SNR(m,1);
        A = diag( [1 A21(n,1)]*SNR(m,1) );
        AA = diag(a);
        
        trial = 0;
        while trial <= TotalTrails
            
            D1 = 2*( rand(P+(K-1)*(P+1),M-P) > 0.5 )-1;
            B1 = [[eye(P);zeros((K-1)*(P+1),P)] D1 ];
            
            if rank(B1)<(M)
                continue
            end
            
            D2 = 2*( rand(P+(K-1)*(P+1),M-P) > 0.5 )-1;
            B2 = [[eye(P);zeros((K-1)*(P+1),P)] D2 ];
            
            if rank(B2)<(M)
                continue
            else
                trial = trial + 1;
            end
            

            b = 2*( rand( P+(P+1)*(K-1), 1 ) > 0.5 )-1;
            rr = SS*AA*b + sigma*randn(P*L,1);

            d1 = pinv(B1)*b;
            d2 = pinv(B2)*b;
            
            %S_new1 = SS*AA*B1;
            S_new1 = SS*AA*B1 + [ zeros(L*P,P) sigma*randn(L*P,M-P)];
            %S_new2 = SS*AA*B2;
            S_new2 = SS*AA*B2 + [ zeros(L*P,P) sigma*randn(L*P,M-P)];
            
            

            %%==============================================================================================
            %%==============================================================================================
            
            d_LS1 = pinv(S_new1'*S_new1)*S_new1'*rr(:,1);
            d_LS2 = pinv(S_new2'*S_new2)*S_new2'*rr(:,1);
            
            A1_LS_LS1 = (D2(1:P,1:M-P)*d_LS2(P+1:M)-D1(1:P,1:M-P)*d_LS1(P+1:M))'*(d_LS1(1:P)-d_LS2(1:P))/norm(D2(1:P,1:M-P)*d_LS2(P+1:M)-D1(1:P,1:M-P)*d_LS1(P+1:M))
            A1_LS_LS2 = norm(d_LS1(1:P)-d_LS2(1:P))/((D2(1:P,1:M-P)*d_LS2(P+1:M)-D1(1:P,1:M-P)*d_LS1(P+1:M))'*(d_LS1(1:P)-d_LS2(1:P)))
            [U_temp,S_temp,V_temp] = svd([(D2(1:P,1:M-P)*d_LS2(P+1:M)-D1(1:P,1:M-P)*d_LS1(P+1:M)) (d_LS1(1:P)-d_LS2(1:P))]);
            eps = min(min(S_temp));
            A1_LS_TLS = (D2(1:P,1:M-P)*d_LS2(P+1:M)-D1(1:P,1:M-P)*d_LS1(P+1:M))'*(d_LS1(1:P)-d_LS2(1:P))/(norm(D2(1:P,1:M-P)*d_LS2(P+1:M)-D1(1:P,1:M-P)*d_LS1(P+1:M))-eps)
            
            b_LS_LS1 = [eye(P)/A1_LS_LS1 B1(1:P,P+1:M)]*d_LS1
            b_LS_LS2 = [eye(P)/A1_LS_LS2 B1(1:P,P+1:M)]*d_LS1
            b_LS_TLS = [eye(P)/A1_LS_TLS B1(1:P,P+1:M)]*d_LS1
            
            dd_LS(trial,1) = norm(d_LS1-d1);            
            %%===============================================================================================%%
            
            [Us, Ss, Vs]=svd([S_new1 rr(:,1)]);
            s_svd=sort(diag(Ss));
            d_TLS1 = pinv( S_new1'*S_new1 - s_svd(1,1)^2*eye((P+1)*K-1) )*S_new1'*rr(:,1);
            
            [Us, Ss, Vs]=svd([S_new2 rr(:,1)]);
            s_svd=sort(diag(Ss));
            d_TLS2 = pinv( S_new2'*S_new2 - s_svd(1,1)^2*eye((P+1)*K-1) )*S_new2'*rr(:,1);

            A1_TLS_LS1 = (D2(1:P,1:M-P)*d_TLS2(P+1:M)-D1(1:P,1:M-P)*d_TLS1(P+1:M))'*(d_TLS1(1:P)-d_TLS2(1:P))/norm(D2(1:P,1:M-P)*d_TLS2(P+1:M)-D1(1:P,1:M-P)*d_TLS1(P+1:M))
            A1_TLS_LS2 = norm(d_TLS1(1:P)-d_TLS2(1:P))/((D2(1:P,1:M-P)*d_TLS2(P+1:M)-D1(1:P,1:M-P)*d_TLS1(P+1:M))'*(d_TLS1(1:P)-d_TLS2(1:P)))
            [U_temp,S_temp,V_temp] = svd([(D2(1:P,1:M-P)*d_TLS2(P+1:M)-D1(1:P,1:M-P)*d_TLS1(P+1:M)) (d_TLS1(1:P)-d_TLS2(1:P))]);
            eps = min(min(S_temp));
            A1_TLS_TLS = (D2(1:P,1:M-P)*d_TLS2(P+1:M)-D1(1:P,1:M-P)*d_TLS1(P+1:M))'*(d_TLS1(1:P)-d_TLS2(1:P))/(norm(D2(1:P,1:M-P)*d_TLS2(P+1:M)-D1(1:P,1:M-P)*d_TLS1(P+1:M))-eps)

            b_TLS_LS1 = [eye(P)/A1_TLS_LS1 B1(1:P,P+1:M)]*d_TLS1
            b_TLS_LS2 = [eye(P)/A1_TLS_LS2 B1(1:P,P+1:M)]*d_TLS1
            b_TLS_TLS = [eye(P)/A1_TLS_TLS B1(1:P,P+1:M)]*d_TLS1
            
            dd_TLS(trial,1) = norm(d_TLS1-d1);            
            %%===============================================================================================%%
            
            temp = Q'*[S_new1 rr(:,1)];
            [Us, Ss, Vs]=svd( temp(P+1:P*L,P+1:M+1) );
            s_svd=sort(diag(Ss));
            d_MLS1 = pinv( S_new1'*S_new1 - diag([zeros(1,P) s_svd(1,1)^2*ones(1, (P+1)*(K-1))]) )*S_new1'*rr(:,1);

            temp = Q'*[S_new2 rr(:,1)];
            [Us, Ss, Vs]=svd( temp(P+1:P*L,P+1:M+1) );
            s_svd=sort(diag(Ss));
            d_MLS2 = pinv( S_new2'*S_new2 - diag([zeros(1,P) s_svd(1,1)^2*ones(1, (P+1)*(K-1))]) )*S_new2'*rr(:,1);
            
            A1_MLS_LS1 = (D2(1:P,1:M-P)*d_MLS2(P+1:M)-D1(1:P,1:M-P)*d_MLS1(P+1:M))'*(d_MLS1(1:P)-d_MLS2(1:P))/norm(D2(1:P,1:M-P)*d_MLS2(P+1:M)-D1(1:P,1:M-P)*d_MLS1(P+1:M))
            A1_MLS_LS2 = norm(d_MLS1(1:P)-d_MLS2(1:P))/((D2(1:P,1:M-P)*d_MLS2(P+1:M)-D1(1:P,1:M-P)*d_MLS1(P+1:M))'*(d_MLS1(1:P)-d_MLS2(1:P)))
            [U_temp,S_temp,V_temp] = svd([(D2(1:P,1:M-P)*d_MLS2(P+1:M)-D1(1:P,1:M-P)*d_MLS1(P+1:M)) (d_MLS1(1:P)-d_MLS2(1:P))]);
            eps = min(min(S_temp));
            A1_MLS_TLS = (D2(1:P,1:M-P)*d_MLS2(P+1:M)-D1(1:P,1:M-P)*d_MLS1(P+1:M))'*(d_MLS1(1:P)-d_MLS2(1:P))/(norm(D2(1:P,1:M-P)*d_MLS2(P+1:M)-D1(1:P,1:M-P)*d_MLS1(P+1:M))-eps)
            
            b_MLS_LS1 = [eye(P)/A1_MLS_LS1 B1(1:P,P+1:M)]*d_MLS1
            b_MLS_LS2 = [eye(P)/A1_MLS_LS2 B1(1:P,P+1:M)]*d_MLS1
            b_MLS_TLS = [eye(P)/A1_MLS_TLS B1(1:P,P+1:M)]*d_MLS1
           
            dd_MLS(trial,1) = norm(d_MLS1-d1);
            %%===============================================================================================%%
            %%=================================================================================================%%

            for p = 1:P
                if ( (b_LS_LS1(p,1)*b(p,1)) < 0 )
                    BER_LS_LS1(m,n) = BER_LS_LS1(m,n) + 1;
                end
                
                if ( (b_TLS_LS1(p,1)*b(p,1)) < 0 )
                    BER_TLS_LS1(m,n) = BER_TLS_LS1(m,n) + 1;
                end
                
                if ( (b_MLS_LS1(p,1)*b(p,1)) < 0 )
                    BER_MLS_LS1(m,n) = BER_MLS_LS1(m,n) + 1;
                end
                
                
                if ( (b_LS_LS2(p,1)*b(p,1)) < 0 )
                    BER_LS_LS2(m,n) = BER_LS_LS2(m,n) + 1;
                end
                
                if ( (b_TLS_LS2(p,1)*b(p,1)) < 0 )
                    BER_TLS_LS2(m,n) = BER_TLS_LS2(m,n) + 1;
                end
                
                if ( (b_MLS_LS2(p,1)*b(p,1)) < 0 )
                    BER_MLS_LS2(m,n) = BER_MLS_LS2(m,n) + 1;
                end
                
                if ( (b_LS_TLS(p,1)*b(p,1)) < 0 )
                    BER_LS_TLS(m,n) = BER_LS_TLS(m,n) + 1;
                end
                
                if ( (b_TLS_TLS(p,1)*b(p,1)) < 0 )
                    BER_TLS_TLS(m,n) = BER_TLS_TLS(m,n) + 1;
                end
                
                if ( (b_MLS_TLS(p,1)*b(p,1)) < 0 )
                    BER_MLS_TLS(m,n) = BER_MLS_TLS(m,n) + 1;
                end
                
                r=rr(p*L-L+1:p*L,1);

                b_MF = sign(w_MF*r);
                if ( (b_MF*b(p,1)) < 0 )
                    BER_MF(m,n) = BER_MF(m,n) + 1;
                end
                
                b_DD = sign(w_DD*r);
                if ( (b_DD*b(p,1)) < 0 )
                    BER_DD(m,n) = BER_DD(m,n) + 1;
                end
                
            end

        end
        
        d_mean_LS(m,n) = mean(dd_LS);
        d_mean_TLS(m,n) = mean(dd_TLS);
        d_mean_MLS(m,n) = mean(dd_MLS);        

        d_std_LS(m,n) = std(dd_LS);
        d_std_TLS(m,n) = std(dd_TLS);
        d_std_MLS(m,n) = std(dd_MLS);        
    end
end

BER_MF = BER_MF/TotalTrails/P
BER_DD = BER_DD/TotalTrails/P

BER_LS_LS1 = BER_LS_LS1/TotalTrails/P
BER_TLS_LS1 = BER_TLS_LS1/TotalTrails/P
BER_MLS_LS1 = BER_MLS_LS1/TotalTrails/P

BER_LS_LS2 = BER_LS_LS2/TotalTrails/P
BER_TLS_LS2 = BER_TLS_LS2/TotalTrails/P
BER_MLS_LS2 = BER_MLS_LS2/TotalTrails/P

BER_LS_TLS = BER_LS_TLS/TotalTrails/P
BER_TLS_TLS = BER_TLS_TLS/TotalTrails/P
BER_MLS_TLS = BER_MLS_TLS/TotalTrails/P

save ABDD20_08172003.mat
pause
%%-------------------------------------------------------------------------------------%
fid = fopen('ASBD15.txt','a');

fprintf(fid,'%f / %f / %f --  %f : %f : %f \n', clock);

fprintf(fid,'d_mean_LS =================================== \n');
fprintf(fid,'%f  %f  %f\n',d_mean_LS);
fprintf(fid,'d_mean_TLS =================================== \n');
fprintf(fid,'%f  %f  %f\n',d_mean_TLS);
fprintf(fid,'d_mean_MLS =================================== \n');
fprintf(fid,'%f  %f %f\n',d_mean_MLS);

fprintf(fid,'d_std_LS =================================== \n');
fprintf(fid,'%f  %f  %f\n',d_std_LS);
fprintf(fid,'d_std_TLS =================================== \n');
fprintf(fid,'%f  %f  %f\n',d_std_TLS);
fprintf(fid,'d_std_MLS =================================== \n');
fprintf(fid,'%f  %f %f\n',d_std_MLS);


fprintf(fid,'BER_MF1 ===================================== \n');
fprintf(fid,'%f  %f %f\n',BER_MF1);
fprintf(fid,'BER_DD1 ===================================== \n');
fprintf(fid,'%f  %f %f\n',BER_DD1);

fprintf(fid,'BER_MF ====================================== \n');
fprintf(fid,'%f  %f %f\n',BER_MF);
fprintf(fid,'BER_DD ====================================== \n');
fprintf(fid,'%f  %f %f\n',BER_DD);

fprintf(fid,'BER_LS ====================================== \n');
fprintf(fid,'%f  %f %f\n',BER_LS);
fprintf(fid,'BER_TLS ===================================== \n');
fprintf(fid,'%f  %f %f\n',BER_TLS);
fprintf(fid,'BER_MLS ===================================== \n');
fprintf(fid,'%f  %f %f\n',BER_MLS);

fclose(fid);
 

%%------------------------------------------------------------------------------------------------%%

SNRlog = log10(SNR).*20;
A21log = log10(A21).*20;

%%------------------------------------------------------------------------------------------------%%

figure(1)
semilogy(SNRlog, d_mean_LS(:,1),'-*',SNRlog, d_mean_TLS(:,1),'-+',SNRlog, d_mean_MLS(:,1),'-x', SNRlog,  d_mean_LS(:,3),'-d', SNRlog,  d_mean_TLS(:,3),'-p', SNRlog,  d_mean_MLS(:,3),'-h')
grid;
title('BER of Different Multiuser Detection Algorithms');
xlabel('SNR (dB)');
ylabel('Mean');
legend('LS with A_2/A_1=0.1','TLS with A_2/A_1=0.1','MLS with A_2/A_1=0.1','LS with A_2/A_1=10','TLS with A_2/A_1=10','MLS with A_2/A_1=10');

figure(2)
semilogy(SNRlog, d_std_LS(:,1),'-*',SNRlog, d_std_TLS(:,1),'-+',SNRlog, d_std_MLS(:,1),'-x', SNRlog,  d_std_LS(:,3),'-d', SNRlog,  d_std_TLS(:,3),'-p', SNRlog,  d_std_MLS(:,3),'-h')
grid;
title('BER of Different Multiuser Detection Algorithms');
xlabel('SNR (dB)');
ylabel('STD');
legend('LS with A_2/A_1=0.1','TLS with A_2/A_1=0.1','MLS with A_2/A_1=0.1','LS with A_2/A_1=10','TLS with A_2/A_1=10','MLS with A_2/A_1=10');

%%-------------------------------------------------------------------------------%%

figure(3)
semilogy(SNRlog, d_mean_LS(:,1),'-*', SNRlog, d_mean_MLS(:,1),'-x', SNRlog,  d_mean_LS(:,3),'-d', SNRlog,  d_mean_MLS(:,3),'-h')
grid;
title('BER of Different Multiuser Detection Algorithms');
xlabel('SNR (dB)');
ylabel('Mean');
legend('LS with A_2/A_1=0.1','MLS with A_2/A_1=0.1','LS with A_2/A_1=10','MLS with A_2/A_1=10');

figure(4)
semilogy(SNRlog, d_std_LS(:,1),'-*', SNRlog, d_std_MLS(:,1),'-x', SNRlog,  d_std_LS(:,3),'-d',  SNRlog,  d_std_MLS(:,3),'-h')
grid;
title('BER of Different Multiuser Detection Algorithms');
xlabel('SNR (dB)');
ylabel('STD');
legend('LS with A_2/A_1=0.1','MLS with A_2/A_1=0.1','LS with A_2/A_1=10','MLS with A_2/A_1=10');

%%---------------------------------------------------------------------------------%%


figure(501)
semilogy(SNRlog, BER_MF(:,1),'-*',SNRlog, BER_DD(:,1),'-+',SNRlog, BER_LS_LS1(:,1),'-x', SNRlog, BER_LS_LS2(:,1),'-o', SNRlog, BER_LS_TLS(:,1),'-d')
ylim([1.0/TotalTrails 0.55])
xlim([min(SNRlog) max(SNRlog)])
grid;
title('BER of Different Multiuser Detection Algorithms');
xlabel('SNR (dB)');
ylabel('Bit-Error-Rate (BER)');
legend('MF with A_2/A_1=10','DD with A_2/A_1=10','LS-LS1 with A_2/A_1=10','LS-LS2 with A_2/A_1=10','LS-TLS with A_2/A_1=10');


figure(503)
semilogy(SNRlog, BER_MF(:,3),'-*',SNRlog, BER_DD(:,3),'-+',SNRlog, BER_LS_LS1(:,3),'-x', SNRlog, BER_LS_LS2(:,3),'-o', SNRlog, BER_LS_TLS(:,3),'-d')
ylim([1.0/TotalTrails 0.55])
xlim([min(SNRlog) max(SNRlog)])
grid;
title('BER of Different Multiuser Detection Algorithms');
xlabel('SNR (dB)');
ylabel('Bit-Error-Rate (BER)');
legend('MF with A_2/A_1=10','DD with A_2/A_1=10','LS-LS1 with A_2/A_1=10','LS-LS2 with A_2/A_1=10','LS-TLS with A_2/A_1=10');


figure(601)
semilogy(SNRlog, BER_MF(:,1),'-*',SNRlog, BER_DD(:,1),'-+',SNRlog, BER_TLS_LS1(:,1),'-x', SNRlog, BER_TLS_LS2(:,1),'-o', SNRlog, BER_TLS_TLS(:,1),'-d')
ylim([1.0/TotalTrails 0.55])
xlim([min(SNRlog) max(SNRlog)])
grid;
title('BER of Different Multiuser Detection Algorithms');
xlabel('SNR (dB)');
ylabel('Bit-Error-Rate (BER)');
legend('MF with A_2/A_1=10','DD with A_2/A_1=10','TLS-LS1 with A_2/A_1=10','TLS-LS2 with A_2/A_1=10','TLS-TLS with A_2/A_1=10');

figure(603)
semilogy(SNRlog, BER_MF(:,3),'-*',SNRlog, BER_DD(:,3),'-+',SNRlog, BER_TLS_LS1(:,3),'-x', SNRlog, BER_TLS_LS2(:,3),'-o', SNRlog, BER_TLS_TLS(:,3),'-d')
ylim([1.0/TotalTrails 0.55])
xlim([min(SNRlog) max(SNRlog)])
grid;
title('BER of Different Multiuser Detection Algorithms');
xlabel('SNR (dB)');
ylabel('Bit-Error-Rate (BER)');
legend('MF with A_2/A_1=10','DD with A_2/A_1=10','TLS-LS1 with A_2/A_1=10','TLS-LS2 with A_2/A_1=10','TLS-TLS with A_2/A_1=10');

figure(701)
semilogy(SNRlog, BER_MF(:,1),'-*',SNRlog, BER_DD(:,1),'-+',SNRlog, BER_MLS_LS1(:,1),'-x', SNRlog, BER_MLS_LS2(:,1),'-o', SNRlog, BER_MLS_TLS(:,1),'-d')
ylim([1.0/TotalTrails 0.55])
xlim([min(SNRlog) max(SNRlog)])
grid;
title('BER of Different Multiuser Detection Algorithms');
xlabel('SNR (dB)');
ylabel('Bit-Error-Rate (BER)');
legend('MF with A_2/A_1=10','DD with A_2/A_1=10','MLS-LS1 with A_2/A_1=10','MLS-LS2 with A_2/A_1=10','MLS-TLS with A_2/A_1=10');

figure(703)
semilogy(SNRlog, BER_MF(:,3),'-*',SNRlog, BER_DD(:,3),'-+',SNRlog, BER_MLS_LS1(:,3),'-x', SNRlog, BER_MLS_LS2(:,3),'-o', SNRlog, BER_MLS_TLS(:,3),'-d')
ylim([1.0/TotalTrails 0.55])
xlim([min(SNRlog) max(SNRlog)])
grid;
title('BER of Different Multiuser Detection Algorithms');
xlabel('SNR (dB)');
ylabel('Bit-Error-Rate (BER)');
legend('MF with A_2/A_1=10','DD with A_2/A_1=10','MLS-LS1 with A_2/A_1=10','MLS-LS2 with A_2/A_1=10','MLS-TLS with A_2/A_1=10');

%%============================================================================================================

figure(8)
semilogy(A21log, BER_MF(1,:),'-*', A21log, BER_DD(1,:),'-o', A21log, BER_LS(1,:),'-d', A21log, BER_TLS(1,:),'-p', A21log, BER_MLS(1,:),'-h')
grid;
title('Near-Far-Resistance');
xlabel('NFR A_2/A_1(dB)');
ylabel('Bit-Error-Rate (BER)');
legend('MF ','DD','LS','TLS','MLS');


figure(11)
semilogy(SNRlog, d_mean_LS(:,1),'-*',SNRlog,  d_mean_LS(:,2),'-s', SNRlog,  d_mean_LS(:,3),'-d')
grid;
title('BER of Different Multiuser Detection Algorithms');
xlabel('SNR (dB)');
ylabel('Mean');
legend('LS with A_2/A_1=0.1','LS with A_2/A_1=1.0','LS with A_2/A_1=10');

figure(12)
semilogy(SNRlog, d_std_LS(:,1),'-*',SNRlog,  d_std_LS(:,2),'-s',SNRlog,  d_std_LS(:,3),'-d')
grid;
title('BER of Different Multiuser Detection Algorithms');
xlabel('SNR (dB)');
ylabel('STD');
legend('LS with A_2/A_1=0.1','LS with A_2/A_1=1.0','LS with A_2/A_1=10');
