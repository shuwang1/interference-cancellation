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
%       1) Start coding: LS, TLS and MLS
%           -- Nov 24, 2004 by Shu Wang, swang@lge.com
%       2) Add LMS, RLS estimator for the amplitude estimation.
%           -- Shu Wang  12/01/2004
%       3) Add seimblind detections. This means the amplitudes are already
%       known for a group of users.
%       
%
%
%
%--------------------------------------------------------------------------

clear all;
close all;

%%---------------------------------------------------
%% Define the system dimensions here.
%%---------------------------------------------------

K0 = 2             %% the number of users
P = 1               %% the detection window. In most cases, especially synchronous cases, it will be set to be 1
L = 64              %% the spreading gain.
LP = L*P
M = K0*4            %% the columns of the blind spreading matrix. 
                    %% Theorectially, the large number, the better performance.  
G = 1               %% the group size. the default is 1. 
IsAsync = 0


%%----------------------------------------------------------
%% define the original spreading sequence vectors/matrix
%%----------------------------------------------------------


%!!! You may feel free to replace this part with other sequences, delay and 
%!!! amplitude definitions.
%% S0 is defined by the original spreading sequences.

%H_walsh = hadamard(L);
%S0 = [H_walsh(1:L/4+3,K0:-1:1); H_walsh(L/4+4:L, 3:K0+2) ]./sqrt(L);                
S0 = ( ( rand(L,K0) > 0.5 ) *2 -1 ) ./sqrt(L);

a0 = ones( K0, 1 ); 
A0 = diag( a0 );      

tao = 0:4:(K0-1)*4;          %% this is the propagation delay vector used 
                             %% in asynchronous cases.

switch  IsAsync 
    %% the spreading sequence matrix for synchronous cases
    case 0        
        K = K0; 
        S = kron( ones(P,1), S0 );           
        a = a0;
        A = diag(a);
        S_G = S(:, 1:G);
    %% initialize the spreading matrix S                  
    case 1
        
        K = 2*K0-1;
        S = zeros( LP, K );            
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
        S_G = S(:, 1:2*G-1);
        
end

R_S = S'*S

R_S_ = inv( R_S );

pause
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

simulationConfig = 1
switch simulationConfig
    case 1
        %
        % Config. setting 1
        % This is to check the near-far resistance
        %
        K = 2
        G = 1
        SNR = 2;
        A21 = [0.1; 0.2; 0.4; 0.6; 0.8; 1; 2; 4; 6; 8; 10];
    case 2
        %
        % Config. setting 2
        % This is to check error rates.
        %
        SNR = [0.1; 0.2; 0.4; 0.6; 0.8; 1; 2; 4; 6; 8;];
        A21 = [0.1; 2; 10];
    case 3
        
end

%%
%% The following decided how to decide the amplitude used for the detection 
%% vector estimation. 
%%

amplitudeFiltering = -1
switch amplitudeFiltering
    case 20
        stepSize = 0.03;
        leakagefactor = 0.5;
        eq_lms = lineareq( M-G, lms(stepSize, leakagefactor) );
        eq_lms.ResetBeforeFiltering = 0;

        eq_LS = eqlms;
        eq_TLS = eqlms;
        eq_MLS = eqlms;
        
    case 21
        stepSize = 0.03;
        bias = 0;
        eq_normlms = lineareq( M-G, normlms(stepSize, bias) );
        eq_normlms.LeakageFactor = 0.5;
        eq_normlms.ResetBeforeFiltering = 0;

        eq_LS = eq_normlms;
        eq_TLS = eq_normlms;
        eq_MLS = eq_normlms;

    case 22
        stepSize = 0.03;
        algtype = 'Sign LMS';
        eq_signlms = lineareq( M-G, lms(stepSize, algtype) );
        eq_signlms.LeakageFactor = 0.5;
        eq_signlms.ResetBeforeFiltering = 0;

        eq_LS = eq_signlms;
        eq_TLS = eq_signlms;
        eq_MLS = eq_signlms;

    case 23
        initstep = 0.005
        incrstep = 0.001
        ministep = 0.005
        maxistep = 0.050
        eq_varlms = lineareq( M-G, varlms(initstep, incrstep, ministep, maxistep ) );
        eq_varlms.ResetBeforeFiltering = 0;

        eq_LS = eq_varlms;
        eq_TLS = eq_varlms;
        eq_MLS = eq_varlms;
        
    case 3
        stepSize = 0.03;
        leakagefactor = 0.5;
        eq_cma = lineareq( M-G, cma(stepSize, leakagefactor) );
        eq_cma.ResetBeforeFiltering = 0;
        
        eq_LS = eq_cma;
        eq_TLS = eq_cma;
        eq_MLS = eq_cma;
        
    case 4
        rlsForgetFactor = 0.5;
        eq_rls = lineareq( M-G, rls(rlsForgetFactor, 0.1) );
        eq_rls.ResetBeforeFiltering = 0;

        eq_LS = eq_rls;
        eq_TLS = eq_rls;
        eq_MLS = eq_rls;
        
end


%%------------------------------------------------------------------------
%%  Play the game.
%%------------------------------------------------------------------------
TotalTrails = 10000

for n = 1:length(A21)
    
    %%
    %% updating ampitudes here
    %%
    
    
    a_ = A21(n,1)*a(G+1:K);
    A = diag([ a(1:G); a_ ])
    
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
        %% ??? When debugging, sigma normally will be set to be zero.
        sigma = 1/SNR(m,1);              

        %%-------------------------
        %% initialization here.
        %%-------------------------
        
        B = [[eye(G); zeros(K-G, G)] ( ( rand(K, M-G+1) > 0.5 )*2 -1 ) ];
        r(:,1:M-G+1) = S*A*B(:,G+1:M+1) + sigma*randn(P*L,M-G+1);

        S_LS = S*A*B(:,1:M);
        G_LS = A(1:G,1:G)*B(1:G,1:M);
        d_LS(1:G,1:M-G+1) = A(1:G,1:G)*B(1:G,G+1:M+1);
        b_LS(1:G,1:M-G+1) = sign( d_LS(1:G,1:M-G+1) );
        a_LS(1:G,1:M-G+1) = abs( d_LS(1:G,1:M-G+1) );
        
        
        S_TLS = S_LS;
        G_TLS = G_LS;
        d_TLS(:,1:M-G+1)  = d_LS(:,1:M-G+1) ;
        b_TLS(:,1:M-G+1) = b_LS(:,1:M-G+1);
        a_TLS(:,1:M-G+1)  = a_LS(:,1:M-G+1);
     
        S_MLS = S_LS;
        G_MLS = G_LS;
        d_MLS(:,1:M-G+1)  = d_LS(:,1:M-G+1) ;
        b_MLS(:,1:M-G+1) = b_LS(:,1:M-G+1);
        a_MLS(:,1:M-G+1)  = a_LS(:,1:M-G+1) ;
        
       
        %%----------------------------------------
        %%  starting simulation trials here.
        %%----------------------------------------

        for t = M-G+1 : TotalTrails
            
            b0 = ( rand(K0,1) > 0.5 ) * 2 -1 ;
            switch ( IsAsync )
                case 0
                    b(:, t) = b0;
                case 1
                    b(2:2:K, t) = b0(2:K0);
                    b(3:2:K, t) = b0(2:K0);
            end
            r(:, t) = S*A*b(:, t) + sigma*randn(P*L,1);
            
            
            %%---------------------------------------------------
            %% Least-Squares Blind MUD
            %%---------------------------------------------------

            %%
            %% updating the so-called blind spreading matrix here.
            %%
            
            S_LS(:,G+1:M) = [ S_LS(:,(G+2):M) r(:, t-1) ];
            switch amplitudeFiltering
                case -1
                    G_LS(:, G+1:M) = [ G_LS(:, (G+2):M) A(1:G,1:G)*b_LS(:, t-1) ];
                case 0
                    G_LS(:, G+1:M) = [ G_LS(:, (G+2):M)  d_LS(:, t-1) ];                    
                case 1
                    G_LS(:, G+1:M) = [ G_LS(:, (G+2):M)  b_LS(:, t-1)*mean(a_LS(:,1:t-1)) ];
                otherwise
                    G_LS(:, G+1:M) = [ G_LS(:, (G+2):M)  abs( equalize( eq_LS, d_LS(:,t-1) ) )*b_LS(:, t-1) ];
            end
            
            %% 
            %% estimate the detection vector using LS algorithm
			%% 
            %% Since S_LS has more rows than columns and may not be of full
            %% rank, then the overdetermined least squares problem may not 
            %% have a unique solution with S_LS of reduced-rank. Two of the
            %% infinitely many solutions are
            %%
            %%      1) f_LS = pinv(S_LS)*r 
            %%      2) f_LS = S_LS\r
            %%
            %% These two are distinguished by the facts that norm(f_LS) is 
            %% the smallest than y has the fewest possible nonzero
            %% components.
            %%
            
            %f_LS = pinv(S_LS)*r(:,t);
            f_LS = S_LS\r(:,t);
            
            
            %% detect the desired users' information bits and estimate its
            %% amplitude.
            d_LS(:, t) = G_LS(1:G,1:M)*f_LS;
            b_LS(:, t) = sign( d_LS(:, t) );
            a_LS(:, t) = abs( d_LS(:, t) );
            
                       
            %%-------------------------------------------------------------
            %% Total Least-Squares Blind MUD
            %%-------------------------------------------------------------

            %%
            %% updating the so-called blind spreading matrix here.
            %%
            
            S_TLS(:, G+1:M) = [ S_TLS(:,(G+2):M) r(:, t-1) ];  
            switch amplitudeFiltering
                case -1
                    G_TLS(:, G+1:M) = [ G_TLS(:, (G+2):M) A(1:G,1:G)*b_TLS(:, t-1) ];
                case 0
                    G_TLS(:, G+1:M) = [ G_TLS(:,(G+2):M) d_TLS(:,t-1) ];
                case 1
                    G_TLS(:, G+1:M) = [ G_TLS(:, (G+2):M)  b_TLS(:, t-1)*mean(a_TLS(:,1:t-1)) ];
                otherwise
                    G_TLS(:, G+1:M) = [ G_TLS(:, (G+2):M)  abs( equalize( eq_TLS, d_TLS(:,t-1) ) )*b_TLS(:, t-1) ];
            end
            
            %% estimate the detection vector using TLS scheme
            %
            % The question here is how to choose the eigenvalue to 
            % compensate the noise, what the rank of [S_TLS r(:,t)] is or
            % maybe what the number of users is.
            %   1) the simplest way is to use the minimum eigvalue
            %   2) use the matlab rank( [S_TLS r(:,t)] )
            %   3) estimate the number of embedded signal sources
            % 

            [U_TLS, D_TLS, V_TLS] = svd([S_TLS r(:,t)],'econ');
            tol = max( size([S_TLS r(:,t)]) )*D_TLS(1,1)*eps;
            %K_TLS = sum( diag(Ds) > tol);            
            K_TLS = M+1;
            
            f_TLS = ( pinv( S_TLS'*S_TLS - (D_TLS(K_TLS, K_TLS)^2)*eye(M) ) * S_TLS' ) * r(:,t);
            
            %% detect the desired users' information bits and estimate its
            %% amplitude.
            d_TLS(:, t) = G_TLS(1:G,1:M) * f_TLS;
            b_TLS(:, t) = sign( d_TLS(:, t) );
            a_TLS(:, t) = abs( d_TLS(:, t) );
            

            %%---------------------------------------------------
            %% Mixed TLS/LS Blind MUD
            %%---------------------------------------------------

            %%
            %% updating the so-called blind spreading matrix here.
            %%
            
            S_MLS(:, G+1:M) = [ S_MLS(:,(G+2):M) r(:, t-1) ];  
            switch amplitudeFiltering
                case -1
                    G_MLS(:, G+1:M) = [ G_MLS(:, (G+2):M) A(1:G,1:G)*b_MLS(:, t-1) ];
                case 0
                    G_MLS(:, G+1:M) = [ G_MLS(:,(G+2):M) d_MLS(:,t-1) ];
                case 1
                    G_MLS(:, G+1:M) = [ G_MLS(:, (G+2):M)  b_MLS(:, t-1)*mean(a_MLS(:,1:t-1)) ];
                otherwise
                    G_MLS(:, G+1:M) = [ G_MLS(:, (G+2):M)  abs( equalize( eq_MLS, d_MLS(:,t-1) ) )*b_MLS(:, t-1) ];
            end
            
            
            %% estimate the detection vector f_MLS using MLS scheme
            [ Qg, Rg ] = qr( S_G );
            S_qr = Qg'*[ S_MLS r(:,t) ];
            [U_MLS, D_MLS, V_MLS] = svd( S_qr(G+1:LP,G+1:M+1), 'econ' );

            tol = max( size( S_qr(G+1:LP,G+1:M+1) ) )*D_MLS(1,1)*eps;
            K_MLS = sum( diag(D_MLS) > tol);            
            % K_MLS = M+1-G;
           
            f_MLS = ( pinv( S_MLS'*S_MLS - (D_MLS(K_MLS,K_MLS)^2)*diag([zeros(1,G) ones(1,M-G)]) ) * S_MLS' ) * r(:,t);            
            
            %% detect the desired users' information bits and estimate its
            %% amplitude.
            d_MLS(:, t) = G_MLS(1:G,1:M) * f_MLS;
            b_MLS(:, t) = sign( d_MLS(:,t) );
            a_MLS(:, t) = abs( d_MLS(:,t) );

            %%-------------------------------------------------------------
            %% multiuser detection using other existing algorithms.
            %%-------------------------------------------------------------
                
            b_MF(:, t) = w_MF' * r(:,t);
            a_MF(:, t) = abs( b_MF(:, t) );
            
            b_DD(:, t) = w_DD' * r(:,t);
            a_DD(:, t) = abs( b_DD(:, t) );

            %%-------------------------------------------------------------
            %% BER counting
            %%-------------------------------------------------------------
            BER_MF(m,n) = BER_MF(m,n) - ( sign(b_MF(1, t))*b(1,t) - 1 )/2;
            BER_DD(m,n) = BER_DD(m,n) - ( sign(b_DD(1, t))*b(1,t) - 1 )/2;

            BER_LS(m,n) = BER_LS(m,n) - ( b_LS(1, t)*b(1,t) - 1 )/2;
            BER_TLS(m,n) = BER_TLS(m,n) - ( b_TLS(1, t)*b(1,t) - 1 )/2;
            BER_MLS(m,n) = BER_MLS(m,n) - ( b_MLS(1, t)*b(1,t) -1 )/2;

        end
        
        %%-----------------------------------------------------------------
        %% Statistical analysis on amplitude estimation.
        %%-----------------------------------------------------------------
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

save(datestr(now, 30))

%%-----------------------------------------------------------------------%%
P1_SNR = SNR.^2;
SNR = log10(SNR).*20;

A2A1 = log10(A21).*20;

%%-------------------------------------------------------------------------

if simulationConfig == 1
    figure(30)
    semilogy(A2A1, BER_MF(1,:),'-*', A2A1, BER_DD(1,:),'-o', A2A1, BER_LS(1,:),'-^', A2A1, BER_TLS(1,:),'-v', A2A1, BER_MLS(1,:),'->','LineWidth',3,'MarkerSize',15)
    grid;
    title('Near-Far Resistance Performance of Different Detectors');
    xlabel('Near-Far Resistance A_2/A_1(dB)','FontSize',14,'FontWeight','bold');
    ylabel('Bit-Error-Rate (BER)','FontSize',14,'FontWeight','bold');
    legend('MF','DD','LS','TLS','MLS');
    return
end

%%-------------------------------------------------------------------------

figure(10001)
plot(1:TotalTrails, a_LS, 1:TotalTrails, a_TLS, 1:TotalTrails, a_MLS)
ylim([0 10]);
grid;
legend('LS','TLS','MLS');

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


%%-------------------------------------------------------------------------


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
%%-------------------------------------------------------------------------


%%=========================================================================
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