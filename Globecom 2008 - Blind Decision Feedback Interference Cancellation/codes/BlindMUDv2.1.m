%
%   
%       Blind Multiuser Detection Framework And Algorithms
%   
%-------------------------------------------------------------------------
%
%   Ideas: 
%       1) propose a blind spreading sequence matrix
%       2) a new blind system model
%       3) and a bunch of algorithms based on this framework.
%
%   Features:
%       1) They are blind.
%       2) They can be one-step or one-shot algorithms. This means that
%       they are very simple and direct.
%       3) They can be adaptive and updatable.
%
%--------------------------------------------------------------------------
%
%   Simulation Purpose:
%       1) verify and show the performance of proposed algorithms, 
%           including LS, TLS, MLS, BLU, MMSE
%       2) compare with MF, DD, MMSE
%       3) 
%
%   major parameters:
%       1) SS -- the blind spreading sequence matrix.
%       2) configSetting -- simulation setting switcher.
%
%--------------------------------------------------------------------------
%
%   History:
%       1) Start coding on Nov 24, 2004 by Shu Wang, swang@lge.com
%       2) Add LMS, RLS estimator for the amplitude estimation.
%           -- Shu Wang  12/01/2004
%       
%
%
%
%---------------------------------------------------------------------------

clear all;
close all;

%%---------------------------------------------------
%% Define the system dimensions here.
%%---------------------------------------------------

K0 = 3                      %% the number of users
P = 1                       %% the detection window. In most cases, especially synchronous cases, it will be set to be 1
L = 16                      %% the spreading gain.
LP = L*P
M = (P+1)*K0-1              %% the columns of the blind spreading matrix. Theorectially, the large number, the better performance.  
G = 1                       %% the group size. the default is 1. 


%%----------------------------------------------------------
%% define the original spreading sequence vectors/matrix
%%----------------------------------------------------------

IsAsync = 0
tao = 0:4:(K0-1)*4;          %% this is the propagation delay vector.

%H_walsh = hadamard(L);
%S0 = [H_walsh(1:L/4+3,K0:-1:1); H_walsh(L/4+4:L, 3:K0+2) ]./sqrt(L);               %% S0 is defined by the original spreading sequences. 
S0 = ( ( rand(L,K0) > 0.5 ) *2 -1 ) ./sqrt(L);

a0 = ones( K0, 1 );                          %!!! You may feel free to replace this part with other sequences and amplitude definitions.
A0 = diag( a0 );                               

switch  IsAsync 
    case 0        
        K = K0; 
        S = kron( ones(P,1), S0 );          %% the spreading sequence matrix for synchronous cases 
        a = a0;
        A = diag(a);
        S_G = S(:, 1:G);
                      
    case 1
        K = 2*K0-1;
        S = zeros( LP, K );            %% initialize the spreading matrix S
        S(:, 1) = kron( ones(P,1), S0(:,1) );
        a = ones(K,1);
        for k = 2:K0

            sk_ = [ S0( L-tao(k)+1:L, k ); zeros(L-tao(k),1) ];
            s_k = [ zeros(tao(k),1); S0( 1:L-tao(k), k ) ];
            a(2*k-2:2*k-1, 1) = [a0(k); a0(k)];
      
            for p = 1:2:P
                S( (p*L-L+1):(p*L), (2*k-2):(2*k-1) ) = [sk_ s_k];
                if p < P
                    S( (p*L+1):(p*L+L), (2*k-2):(2*k-1) ) = [s_k sk_];
                end
            end
            
        end
        A = diag( a );
        S_G = S(:, 1:G);
        
end

R_S = S'*S
R_S_ = inv( R_S );

%%--------------------------------------------------------
%% Construct MF receiver and several conventional MUDs.
%%--------------------------------------------------------

w_MF  = S(:,1);
w_DD  = ( R_S_ (1,:)*S' )';


%%-------------------------------------------------------
%% Define simulation environment parameters AND start
%% the simulations.
%%-------------------------------------------------------

sigma = 1;

simulationConfig = 2
switch simulationConfig
    case 1
        %
        % Config. setting 1
        % This is to check the near-far resistance
        %
        SNR = 2;
        A21 = [0.1; 0.2; 0.4; 0.6; 0.8; 1; 2; 4; 6; 8; 10];
    case 2
        %
        % Config. setting 2
        % This is to check error rates.
        %
        SNR = [0.1; 0.2; 0.4; 0.6; 0.8; 1; 2; 3; 4; 5];
        A21 = [0.1; 1; 10];
end

rlsForgetFactor = 0.5;

eqrls_LS = lineareq( M, rls(rlsForgetFactor, 0.1) );
eqrls_LS.ResetBeforeFiltering = 0;

eq_current = eqrls_LS;
%%-------------------------------------------------------
%%  Play the game.
%%-------------------------------------------------------
TotalTrails = 1000

for n = 1:length(A21)
    
    %%
    %% updating ampitudes here
    %%
    
    
    a1(2:K,1) = A21(n,1)*a(2:K);
    A = diag([a(1); a1(2:K)])
    
    for m = 1:length(SNR)
        
        %
        % initalizing some statistical variables
        %
        BER_MF(m,n) = 0;
        BER_DD(m,n) = 0;
        BER_MMSE(m,n) = 0;
        
        BER_LS(m,n) = 0;
        BER_TLS(m,n) = 0;
        BER_MLS(m,n) = 0;
        BER_BLU(m,n) = 0;
        BER_BMMS(m,n) = 0;
        
        %
        %
        %
        A_MF_mean(m,n) = 0;
        A_DD_mean(m,n) = 0;
        A_MMSE_mean(m,n) = 0;
        
        A_LS_mean(m,n) = 0;
        A_TLS_mean(m,n) = 0;
        A_MLS_mean(m,n) = 0;
        A_BLU_mean(m,n) = 0;
        A_BMMS_mean(m,n) = 0;
        
        A_MF_std(m,n) = 0;
        A_DD_std(m,n) = 0;
        A_MMSE_std(m,n) = 0;
        
        A_LS_std(m,n) = 0;
        A_TLS_std(m,n) = 0;
        A_MLS_std(m,n) = 0;
        A_BLU_std(m,n) = 0;
        A_BMMS_std(m,n) = 0;
        
        %%
        %% updating SNR here
        %%
        sigma = 1/SNR(m,1);             %% ??? When debugging, sigma normally will be set to be zero. 

        %%-------------------------
        %% initialization here.
        %%-------------------------
        
        B = A*( ( rand(K, M) > 0.5 )*2 -1 );
        S_LS = S*B;
        G_LS = B(1:G,1:M);
        a_LS = a(1:G);
        
        S_TLS = S_LS;
        G_TLS = G_LS;
        a_TLS = a_LS;
        
        S_MLS = S_LS;
        G_MLS = G_LS;
        a_MLS = a_LS;
        
        b0 = ( rand(K0,1) > 0.5 ) * 2 -1;        
        switch ( IsAsync )            
            case 0
                b = b0;
            case 1
                b(2:2:K) = b0(2:K0);
                b(3:2:K) = b0(2:K0);
        end
        r = S*A*b + sigma*randn(P*L,1);
        d_LS(:,1) = A(1:G,1:G)*b(1:G);
        d_TLS(:,1) = A(1:G,1:G)*b(1:G);
        d_MLS(:,1) = A(1:G,1:G)*b(1:G);
               
        
        %%----------------------------------------
        %%  starting simulation trials here.
        %%----------------------------------------

        trial = 2;
        
        
        while trial <= TotalTrails
            
            r_ = r;

            b0 = ( rand(K0,1) > 0.5 ) * 2 -1;
            switch ( IsAsync )
                case 0
                    b = b0;
                case 1
                    b(2:2:K) = b0(2:K0);
                    b(3:2:K) = b0(2:K0);
            end
            r = S*A*b + sigma*randn(P*L,1);
            
            
            %%---------------------------------------------------
            %% Least-Squares Blind MUD
            %%---------------------------------------------------

            %%
            %% updating the so-called blind spreading matrix here.
            %%
            
            S_LS = [ S_G S_LS(:,(G+2):M) r_ ];
            G_LS = [ eye(G) G_LS(:, (G+2):M)  d_LS(:,trial-1) ];
           
            f_LS = ( pinv(S_LS'*S_LS) * S_LS' ) * r;             %% estimate the detection vector using LS algorithm
            d_LS(:, trial) = equalize( eq_current, G_LS*f_LS );
            b_LS = sign( d_LS(:, trial) );
            a_LS(:, trial) = abs( d_LS(:, trial) );
            
                       
            %%---------------------------------------------------
            %% Total Least-Squares Blind MUD
            %%---------------------------------------------------

            %%
            %% updating the so-called blind spreading matrix here.
            %%
            
            S_TLS = [ S_G S_TLS(:,(G+2):M) r_ ];  
            G_TLS = [ eye(G) G_TLS(:, (G+2):M)  sign(d_TLS(:,trial-1))*mean(a_TLS(1:G,:)) ];
            
            [Us, Ds, Vs] = svd([S_TLS r]);
            d_svd = sort( diag(Ds) );
            
            f_TLS = ( pinv( S_TLS'*S_TLS - (d_svd(1,1)^2)*eye(M) ) * S_TLS' ) * r;            %% estimate the detection vector using TLS scheme
            d_TLS(:, trial) = G_TLS * f_TLS;
            b_TLS = sign( d_TLS(:, trial) );
            a_TLS(:, trial) = abs( d_TLS(:, trial) );
            

            %%---------------------------------------------------
            %% Mixed TLS/LS Blind MUD
            %%---------------------------------------------------

            %%
            %% updating the so-called blind spreading matrix here.
            %%
            
            S_MLS = [ S_G S_MLS(:,(G+2):M) r_ ];  
            G_MLS = [ eye(G) G_MLS(:, (G+2):M)  sign(d_MLS(:,trial-1))*mean(a_MLS(1:G,:)) ];
            
            [ Qg, Rg ] = qr( S_G );
            S_qr = Qg'*[ S_MLS r ];
            [Ug, Dg, Vg] = svd( S_qr(G+1:LP,G+1:M+1) );
            d_svd = sort( diag(Dg) );
            
            f_MLS = ( pinv( S_MLS'*S_MLS - (d_svd(1,1)^2)*diag([0 ones(1,M-G)]) ) * S_MLS' ) * r;            %% estimate the detection vector using MLS scheme
            d_MLS(:, trial) = G_MLS * f_MLS;
            b_MLS = sign( d_MLS(:,trial) );
            a_MLS(:, trial) = abs( d_MLS(:,trial) );

            %=================================================%
                
            b_MF = w_MF' * r(:,1);
            a_MF(1, trial) = abs( b_MF );
            
            b_DD = w_DD' * r(:,1);
            a_DD(1, trial) = abs( b_DD );
            
            BER_MF(m,n) = BER_MF(m,n) - ( sign(b_MF(1,1))*b(1,1) - 1 )/2;
            BER_DD(m,n) = BER_DD(m,n) - ( sign(b_DD(1,1))*b(1,1) - 1 )/2;

            BER_LS(m,n) = BER_LS(m,n) - ( b_LS(1,1)*b(1,1) - 1 )/2;
            BER_TLS(m,n) = BER_TLS(m,n) - ( b_TLS(1,1)*b(1,1) - 1 )/2;
            BER_MLS(m,n) = BER_MLS(m,n) - ( b_MLS(1,1)*b(1,1) -1 )/2;

            trial = trial + 1;
        end
        
        A_LS_mean(m,n) = mean(a_LS(1,:));
        A_TLS_mean(m,n) = mean(a_TLS(1,:));
        A_MLS_mean(m,n) = mean(a_MLS(1,:));
        A_DD_mean(m,n) = mean(a_DD(1,:));
        A_MF_mean(m,n) = mean(a_MF(1,:));        
        
        A_LS_std(m,n) = std(a_LS(1,:));
        A_TLS_std(m,n) = std(a_TLS(1,:));
        A_MLS_std(m,n) = std(a_MLS(1,:));        
        A_DD_std(m,n) = std(a_DD(1,:));
        A_MF_std(m,n) = std(a_MF(1,:));        

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

%%-------------------------------------------------------------------------
%%-------------------------------------------------------------------------

figure(5)
semilogy(SNR, A_MF_std(:,1),'-s', SNR, A_DD_std(:,1),'-o',SNR, A_LS_std(:,1),'-^', SNR, A_TLS_std(:,1),'-v', SNR , A_MLS_std(:,1),'-<' ,'LineWidth',2,'MarkerSize',5)
grid;
title('The Perofmance of The Amplitude Estimation Schemes with A_2/A_1=0.1');
xlabel('SNR (dB)','FontSize',14,'FontWeight','bold');
ylabel('STD of the estimation of A_1','FontSize',14,'FontWeight','bold');
legend('MF','DD','LS','TLS','MLS');

figure(6)
semilogy(SNR, A_MF_std(:,3),'-s', SNR, A_DD_std(:,3),'-o',SNR, A_LS_std(:,3),'-^', SNR, A_TLS_std(:,3),'-v', SNR , A_MLS_std(:,3),'-<','LineWidth',2,'MarkerSize',5)
grid;
title('The Perofmance of The Amplitude Estimation Schemes with A_2/A_1=10');
xlabel('SNR (dB)','FontSize',14,'FontWeight','bold');
ylabel('STD of the estimation of A_1','FontSize',14,'FontWeight','bold');
legend('MF','DD','LS','TLS','MLS');

figure(7)
semilogy(SNR, A_MF_mean(:,1),'-s', SNR, A_DD_mean(:,1),'-o',SNR, A_LS_mean(:,1),'-^', SNR, A_TLS_mean(:,1),'-v', SNR , A_MLS_mean(:,1),'-<','LineWidth',2,'MarkerSize',5)
grid;
title('The Perofmance of The Amplitude Estimation Schemes with A_2/A_1=0.1');
xlabel('SNR (dB)','FontSize',14,'FontWeight','bold');
ylabel('Mean of the estimation of A_1','FontSize',14,'FontWeight','bold');
legend('MF','DD','LS','TLS','MLS');

figure(8)
semilogy(SNR, A_MF_mean(:,3),'-s', SNR, A_DD_mean(:,3),'-o',SNR, A_LS_mean(:,3),'-^', SNR, A_TLS_mean(:,3),'-v', SNR , A_MLS_mean(:,3),'-<','LineWidth',2,'MarkerSize',5)
grid;
title('The Perofmance of The Amplitude Estimation Schemes with A_2/A_1=10');
xlabel('SNR (dB)','FontSize',14,'FontWeight','bold');
ylabel('Mean of the estimation of A_1','FontSize',14,'FontWeight','bold');
legend('MF','DD','LS','TLS','MLS');


figure(13)
semilogy(SNR, A_MF_std(:,1),'-+', SNR, A_MF_std(:,3),'-*', SNR, A_DD_std(:,1),'-o', SNR, A_DD_std(:,3),'-s',SNR, A_LS_std(:,1),'-^',SNR, A_LS_std(:,3),'-v',SNR, A_MLS_std(:,1),'-<',SNR, A_MLS_std(:,3),'->','LineWidth',2,'MarkerSize',5)
grid;
title('The Perofmance of The Amplitude Estimation Schemes with A_2/A_1=0.1');
xlabel('SNR (dB)','FontSize',14,'FontWeight','bold');
ylabel('STD of the estimation of A_1','FontSize',14,'FontWeight','bold');
legend('MF with A_2/A_1=0.1','MF with A_2/A_1=10','DD with A_2/A_1=0.1','DD with A_2/A_1=10','LS with A_2/A_1=0.1','LS with A_2/A_1=10','MLS with A_2/A_1=0.1','MLS with A_2/A_1=10');

figure(14)
semilogy(SNR,ones(size(SNR)),'-.',SNR, A_MF_mean(:,1),'-+',SNR, A_MF_mean(:,3),'-*', SNR, A_DD_mean(:,1),'-o', SNR, A_DD_mean(:,3),'-s',SNR, A_LS_mean(:,1),'-^',SNR, A_LS_mean(:,3),'-v',SNR, A_MLS_mean(:,1),'-<',SNR, A_MLS_mean(:,3),'->','LineWidth',2,'MarkerSize',5)
grid;
title('The Perofmance of The Amplitude Estimation Schemes with A_2/A_1=0.1');
xlabel('SNR (dB)','FontSize',14,'FontWeight','bold');
ylabel('Mean of the estimation of A_1','FontSize',14,'FontWeight','bold');
legend('The orginal power','MF with A_2/A_1=0.1','MF with A_2/A_1=10','DD with A_2/A_1=0.1','DD with A_2/A_1=10','LS with A_2/A_1=0.1','LS with A_2/A_1=10','MLS with A_2/A_1=0.1','MLS with A_2/A_1=10');


figure(15)
semilogy(SNR, A_MF_std(:,1),'-+', SNR, A_MF_std(:,2),'-*', SNR, A_MF_std(:,3),'-x', SNR, A_DD_std(:,1),'-d', SNR, A_DD_std(:,2),'-o', SNR, A_DD_std(:,3),'-s',SNR, A_LS_std(:,1),'-^',SNR, A_LS_std(:,2),'-v',SNR, A_LS_std(:,3),'-<','LineWidth',2,'MarkerSize',5)
grid;
title('The Perofmance of The Amplitude Estimation Schemes with A_2/A_1=0.1');
xlabel('SNR (dB)','FontSize',14,'FontWeight','bold');
ylabel('STD of the estimation of A_1','FontSize',14,'FontWeight','bold');
legend('MF with A_2/A_1=0.1','MF with A_2/A_1=1.0','MF with A_2/A_1=10','DD with A_2/A_1=0.1','DD with A_2/A_1=1.0','DD with A_2/A_1=10','LS with A_2/A_1=0.1','LS with A_2/A_1=1.0','LS with A_2/A_1=10');

figure(16)
semilogy(SNR,ones(size(SNR)),'-.',SNR, A_MF_mean(:,1),'-+',SNR, A_MF_mean(:,2),'-x',SNR, A_MF_mean(:,3),'-*',SNR, A_DD_mean(:,1),'-d', SNR, A_DD_mean(:,2),'-o', SNR, A_DD_mean(:,3),'-s',SNR, A_LS_mean(:,1),'-^',SNR, A_LS_mean(:,2),'-v',SNR, A_LS_mean(:,3),'-<','LineWidth',2,'MarkerSize',5)
grid;
title('The Perofmance of The Amplitude Estimation Schemes with A_2/A_1=0.1');
xlabel('SNR (dB)','FontSize',14,'FontWeight','bold');
ylabel('Mean of the estimation of A_1','FontSize',14,'FontWeight','bold');
legend('The orginal power','MF with A_2/A_1=0.1','MF with A_2/A_1=1.0','MF with A_2/A_1=10','DD with A_2/A_1=0.1','DD with A_2/A_1=1.0','DD with A_2/A_1=10','LS with A_2/A_1=0.1','LS with A_2/A_1=1.0','LS with A_2/A_1=10');


%%------------------------------------------------------------------------------------

figure(27)
semilogy(SNR, BER_MF(:,1),'-*',SNR, BER_MF(:,2),'-+',SNR, BER_MF(:,3),'-x', SNR, BER_DD(:,2),'-o', SNR, BER_LS(:,1),'-d', SNR, BER_LS(:,2),'-p', SNR, BER_LS(:,3),'-h','LineWidth',2,'MarkerSize',5)
ylim([1.0/TotalTrails 0.55])
xlim([min(SNR) max(SNR)])
grid;
title('BER of Different Multiuser Detection Algorithms');
xlabel('SNR (dB)','FontSize',14,'FontWeight','bold');
ylabel('Bit-Error-Rate BER','FontSize',14,'FontWeight','bold');
legend('MF with A_2/A_1=0.1','MF with A_2/A_1=1.0','MF with A_2/A_1=10','DD with A_2/A_1=1.0','LS with A_2/A_1=0.1','LS with A_2/A_1=1.0','LS with A_2/A_1=10');

figure(28)
semilogy(SNR, BER_MF(:,1),'-*',SNR, BER_MF(:,2),'-+',SNR, BER_MF(:,3),'-x', SNR, BER_DD(:,2),'-o', SNR, BER_TLS(:,1),'-d', SNR, BER_TLS(:,2),'-p', SNR, BER_TLS(:,3),'-h','LineWidth',2,'MarkerSize',5)
ylim([1.0/TotalTrails 0.55])
xlim([min(SNR) max(SNR)])
grid;
title('BER of Different Multiuser Detection Algorithms');
xlabel('SNR (dB)','FontSize',14,'FontWeight','bold');
ylabel('Bit-Error-Rate BER','FontSize',14,'FontWeight','bold');
legend('MF with A_2/A_1=0.1','MF with A_2/A_1=1.0','MF with A_2/A_1=10','DD with A_2/A_1=1.0','TLS with A_2/A_1=0.1','TLS with A_2/A_1=1.0','TLS with A_2/A_1=10');

figure(29)
semilogy(SNR, BER_MF(:,1),'-*',SNR, BER_MF(:,2),'-+',SNR, BER_MF(:,3),'-x', SNR, BER_DD(:,2),'-o', SNR, BER_MLS(:,1),'-d', SNR, BER_MLS(:,2),'-p', SNR, BER_MLS(:,3),'-h','LineWidth',2,'MarkerSize',5)
ylim([1.0/TotalTrails 0.55])
xlim([min(SNR) max(SNR)])
grid;
title('BER of Different Multiuser Detection Algorithms');
xlabel('SNR (dB)','FontSize',14,'FontWeight','bold');
ylabel('Bit-Error-Rate BER','FontSize',14,'FontWeight','bold');
legend('MF with A_2/A_1=0.1','MF with A_2/A_1=1.0','MF with A_2/A_1=10','DD with A_2/A_1=1.0','MLS with A_2/A_1=0.1','MLS with A_2/A_1=1.0','MLS with A_2/A_1=10');

return
%%------------------------------------------------------------------------------------------
figure(30)
semilogy(A2A1, BER_MF(1,:),'-*', A2A1, BER_DD(1,:),'-o', A2A1, BER_LS(1,:),'-^', A2A1, BER_TLS(1,:),'-v','LineWidth',2,'MarkerSize',5)
grid;
title('Near-Far Resistance Performance of Different Detectors');
xlabel('Near-Far Resistance A_2/A_1(dB)','FontSize',14,'FontWeight','bold');
ylabel('Bit-Error-Rate (BER)','FontSize',14,'FontWeight','bold');
legend('MF','DD','LS','TLS');

figure(31)
semilogy(A2A1, BER_MF(1,:),'-*', A2A1, BER_DD(1,:),'-o', A2A1, BER_LS(1,:),'-^', A2A1, BER_MLS(1,:),'-v','LineWidth',2,'MarkerSize',5)
grid;
title('Near-Far Resistance Performance of Different Detectors');
xlabel('Near-Far Resistance A_2/A_1(dB)','FontSize',14,'FontWeight','bold');
ylabel('Bit-Error-Rate (BER)','FontSize',14,'FontWeight','bold');
legend('MF','DD','LS','MLS');

figure(32)
semilogy(A2A1, BER_MF(1,:),'-*', A2A1, BER_DD(1,:),'-o', A2A1, BER_LS(1,:),'-^','LineWidth',2,'MarkerSize',5)
grid;
title('Near-Far Resistance Performance of Different Detectors');
xlabel('Near-Far Resistance A_2/A_1(dB)','FontSize',14,'FontWeight','bold');
ylabel('Bit-Error-Rate (BER)','FontSize',14,'FontWeight','bold');
legend('MF','DD','LS');

%%==================================================================================================================%
figure(43)
semilogy(A2A1,ones(size(A2A1)),'-.',A2A1, A_MF_mean(1,:),'-*', A2A1, A_DD_mean(1,:),'-o', A2A1, A_LS_mean(1,:),'-^', A2A1, A_MLS_mean(1,:),'-v','LineWidth',2,'MarkerSize',5)
grid;
title('Near-Far Resistance Performance of Different Detectors');
xlabel('Near-Far Resistance A_2/A_1(dB)','FontSize',14,'FontWeight','bold');
ylabel('Mean of the estimation of A_1','FontSize',14,'FontWeight','bold');
legend('The original power','MF','DD','LS','MLS');

figure(44)
semilogy(A2A1, A_MF_std(1,:),'-*', A2A1, A_DD_std(1,:),'-o', A2A1, A_LS_std(1,:),'-^', A2A1, A_MLS_std(1,:),'-v','LineWidth',2,'MarkerSize',5)
grid;
title('Near-Far Resistance Performance of Different Detectors');
xlabel('Near-Far Resistance A_2/A_1(dB)','FontSize',14,'FontWeight','bold');
ylabel('STD of the estimation of A_1','FontSize',14,'FontWeight','bold');
legend('MF','DD','LS','MLS');


figure(45)
semilogy(A2A1,ones(size(A2A1)),'-.',A2A1, A_MF_mean(1,:),'-*', A2A1, A_DD_mean(1,:),'-o', A2A1, A_LS_mean(1,:),'-^','LineWidth',2,'MarkerSize',5)
grid;
title('Near-Far Resistance Performance of Different Detectors');
xlabel('Near-Far Resistance A_2/A_1(dB)','FontSize',14,'FontWeight','bold');
ylabel('Mean of the estimation of A_1','FontSize',14,'FontWeight','bold');
legend('The original power','MF','DD','LS');

figure(46)
semilogy(A2A1, A_MF_std(1,:),'-*', A2A1, A_DD_std(1,:),'-o', A2A1, A_LS_std(1,:),'-^','LineWidth',2,'MarkerSize',5)
grid;
title('Near-Far Resistance Performance of Different Detectors');
xlabel('Near-Far Resistance A_2/A_1(dB)','FontSize',14,'FontWeight','bold');
ylabel('STD of the estimation of A_1','FontSize',14,'FontWeight','bold');
legend('MF','DD','LS');