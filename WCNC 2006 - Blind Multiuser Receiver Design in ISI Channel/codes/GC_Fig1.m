BER_MF =[    0.2236    0.2286    0.3826;
    0.1758    0.1802    0.3682;
    0.1284    0.1464    0.3856;
    0.1012    0.1084    0.3798;
    0.0574    0.0792    0.3666;
    0.0316    0.0516    0.3704;
    0.0166    0.0404    0.3800;
    0.0060    0.0292    0.3676;
    0.0022    0.0246    0.3604;
    0.0002    0.0198    0.3712;
         0    0.0186    0.3650;
         0    0.0146    0.3660]
         
         
         
BER_DD =[  0.2374    0.2426    0.2362
    0.1938    0.1848    0.2012;
    0.1244    0.1378    0.1406;
    0.0934    0.0842    0.0854;
    0.0394    0.0446    0.0392;
    0.0140    0.0134    0.0134;
    0.0032    0.0044    0.0028;
    0.0004    0.0006         0;
    0.0002         0         0;
         0         0         0;
         0         0         0;
         0         0         0 ]        
         
         
BER_LS = [ 0.2308    0.2334    0.3846;
    0.1884    0.1914    0.3700;
    0.1190    0.1428    0.3758;
    0.0820    0.0928    0.3610;
    0.0356    0.0552    0.3492;
    0.0126    0.0280    0.3386;
    0.0020    0.0128    0.3162;
    0.0004    0.0050    0.2882;
         0    0.0012    0.2504;
         0    0.0004    0.2188;
         0    0.0002    0.1096;
         0         0    0.2154]
         
         
BER_TLS = [ 0.5018    0.4962    0.5024;
    0.4936    0.4934    0.5000;
    0.4988    0.5016       0.5000;
    0.4964    0.4986       0.5000;
    0.4850    0.5038       0.5000;
    0.5068    0.4880       0.5000;
    0.4930    0.5090       0.5000;
    0.5070    0.5048    0.5074;
    0.4932    0.5070    0.4990;
    0.4852    0.5024       0.5000;
    0.4926    0.4990    0.4936;
    0.5016    0.5030       0.5000 ]        
    
    
BER_MLS = [  0.4846    0.4900    0.4934;
    0.4756    0.5102    0.5134;
    0.5008    0.5022    0.4880;
    0.5000    0.4906    0.5012;
    0.4904    0.5012    0.4892;
    0.4996    0.4944    0.4882;
    0.5002    0.5002    0.4896;
    0.5062    0.5038    0.4946;
    0.5030    0.4856    0.4942;
    0.4972    0.4956    0.5050;
    0.4888    0.5022    0.4948;
    0.4932    0.5018    0.4906 ]   
    
    
    
BER_MF1 = [ 2.1735e-001  2.2630e-001  3.8575e-001;
  1.6820e-001  1.7610e-001  3.7840e-001;
  1.1705e-001  1.3465e-001  3.8570e-001;
  7.5500e-002  9.7650e-002  3.8525e-001;
  3.9600e-002  6.2600e-002  3.8055e-001;
  1.5400e-002  4.2050e-002  3.7985e-001;
  5.4500e-003  2.6100e-002  3.8625e-001;
  1.3500e-003  1.6950e-002  3.8040e-001;
  1.0000e-004  1.2450e-002  3.7955e-001;
            0  8.8000e-003  3.7995e-001;
            0  7.6500e-003  3.7775e-001;
            0  5.7500e-003  3.7085e-001 ]        
         
BER_DD1 = [ 2.4215e-001  2.4150e-001  2.4025e-001;
  1.8660e-001  1.8795e-001  1.8245e-001;
  1.2785e-001  1.3270e-001  1.3150e-001;
  8.0750e-002  8.4050e-002  7.8500e-002;
  3.6650e-002  3.7650e-002  4.1450e-002;
  1.2250e-002  1.2850e-002  1.1550e-002;
  1.9500e-003  2.7500e-003  2.5500e-003;
  4.0000e-004  2.0000e-004  1.0000e-004;
            0            0            0;
            0            0            0;
            0            0            0;
            0            0            0 ]        
            
BER_LS1 = [  2.2535e-001  2.3095e-001  3.3585e-001 ;
  1.6955e-001  1.8080e-001  2.9380e-001;
  1.1070e-001  1.3715e-001  2.4890e-001;
  6.5750e-002  9.3500e-002  1.8035e-001;
  2.7000e-002  5.7100e-002  1.1095e-001;
  7.4500e-003  3.3350e-002  5.0150e-002;
  1.0500e-003  1.8500e-002  1.3950e-002;
  1.0000e-004  9.2500e-003  2.3000e-003;
            0  3.5000e-003            0;
            0  1.0000e-003            0;
            0  5.0000e-005            0;
            0            0            0 ]           
            
TotalTrails = 5000             

        SNR = 10.^(-0.1:0.1:1)';

figure(100)
        
        
        subplot(1,2,1);
        semilogy(SNR, BER_MF(:,2),'-+',SNR, BER_MF(:,3),'-x', SNR, BER_DD(:,2),'-o', SNR, BER_LS(:,2),'-<', SNR, BER_LS(:,3), '->',SNR, BER_TLS(:,2),'-^', SNR, BER_TLS(:,3),'-v', SNR, BER_MLS(:,2),'-s', SNR, BER_MLS(:,3),'-d','LineWidth',4,'MarkerSize',10)
        ylim([1.0/TotalTrails 0.55])
        xlim([min(SNR) max(SNR)])
        grid;
        title('BER of Different Multiuser Detection Algorithms');
        xlabel('SNR (dB)','FontSize',24,'FontWeight','bold');
        ylabel('BER','FontSize',24,'FontWeight','bold');
        
legend(...
  {'MF with A_2/A_1=1.0','MF with A_2/A_1=10','DD with A_2/A_1=1.0','LS with A_2/A_1=1.0','LS with A_2/A_1=10','TLS with A_2/A_1=1.0','TLS with A_2/A_1=10','MLS with A_2/A_1=1.0','MLS with A_2/A_1=10'},...
  'FontSize',14,...
  'FontWeight','bold',...
  'Location','Best');
          
        
TotalTrails = 20000        
        subplot(1,2,2);
        
semilogy(SNR, BER_MF1(:,1),'-*',SNR, BER_MF1(:,2),'-+',SNR, BER_MF1(:,3),'-x', SNR, BER_DD1(:,2),'-o', SNR, BER_LS1(:,1),'-s', SNR, BER_LS1(:,2),'-<', SNR, BER_LS1(:,3),'->','LineWidth',4,'MarkerSize',10)
ylim([1.0/TotalTrails 0.55])
xlim([min(SNR) max(SNR)])
grid;
title('BER of Different Multiuser Detection Algorithms');
xlabel('SNR (dB)','FontSize',24,'FontWeight','bold');
ylabel('BER','FontSize',24,'FontWeight','bold');
legend(...
  {'MF with A_2/A_1=0.1','MF with A_2/A_1=1.0','MF with A_2/A_1=10','DD with A_2/A_1=1.0','LS with A_2/A_1=0.1','LS with A_2/A_1=1.0','LS with A_2/A_1=10'},...
  'FontSize',14,...
  'FontWeight','bold',...
  'Location','Best');
  
        
            