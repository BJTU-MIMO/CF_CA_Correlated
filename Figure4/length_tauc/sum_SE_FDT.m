
clear;

%parpool('local',16);

%% Define simulation setup

nbrOfSetups =200; %200;

Bandwith=20e6;

L = 100;

K =20;

N=2;

tau_p =10;

aa=[0.01 0.05];

train=tau_p:-1:1;

P_pilot = 100;
P_u=100;
P_d=200;

ASD=30;

FDT=0:0.0005:0.004;

tau_c2=floor(2.4048./(2*pi*0.004));
tau_d2=tau_c2-tau_p;
data_u2=0:1:(tau_d2-1);
data_d2=0:1:(tau_d2-1);

average_sumSE=zeros(length(FDT),length(aa),nbrOfSetups);
average_sumSE2=zeros(length(FDT),nbrOfSetups);

%%
tic;
for n = 1:nbrOfSetups    
     disp(['Setup------------------- ' num2str(n) ' out of ' num2str(nbrOfSetups)]); 
     
             [gainOverNoisedB,R,pilotIndex] = functionSetup(K,L,N,tau_p,ASD,Bandwith);
         for s=1:length(FDT)
             disp(['FDT------ ' num2str(s) ' out of ' num2str(length(FDT))]);
             
             rhotaup=besselj(0,2*pi*FDT(s)*train);
             [gamma_kl,gamma_kil,tracQ] = functionChannelEstimates(R,K,L,N,P_pilot,pilotIndex,rhotaup);
             
             % full power
             mu1=1./repmat(sum(tracQ(:,:),1),K,1);  
             
             for r=1:length(aa)
                          
              tau_c =tau_p/aa(r);
              tau_d=tau_c-tau_p;
              data_u=0:1:(tau_d-1);
              data_d=0:1:(tau_d-1);
             rhon_u=besselj(0,2*pi*FDT(s)*data_u); 
             rhon_d=besselj(0,2*pi*FDT(s)*data_d); 
            
             [SE_MR_analytical_LSFD] = functionComputeSE_AP_uplink_analytical_LSFD(R,gamma_kl,gamma_kil,K,L,tau_c,tau_d,P_u,pilotIndex,rhon_u);          
             [SE_MR_analytical_coherent] = functionComputeSE_AP_downlink_analytical_coherent(R,gamma_kl,gamma_kil,K,L,tau_c,tau_d,P_d,mu1,pilotIndex,rhon_d);             
             [sum_SE] = functionCompute_sum_SE(SE_MR_analytical_LSFD,SE_MR_analytical_coherent);
             average_sumSE(s,r,n)=sum_SE;
             
             end
             
             rhon_u2=besselj(0,2*pi*FDT(s)*data_u2); 
             rhon_d2=besselj(0,2*pi*FDT(s)*data_d2); 
             
             [SE_MR_analytical_LSFD] = functionComputeSE_AP_uplink_analytical_LSFD(R,gamma_kl,gamma_kil,K,L,tau_c2,tau_d2,P_u,pilotIndex,rhon_u2);          
             [SE_MR_analytical_coherent] = functionComputeSE_AP_downlink_analytical_coherent(R,gamma_kl,gamma_kil,K,L,tau_c2,tau_d2,P_d,mu1,pilotIndex,rhon_d2);             
             [sum_SE] = functionCompute_sum_SE(SE_MR_analytical_LSFD,SE_MR_analytical_coherent);
             average_sumSE2(s,n)=sum_SE;
             
                  
         end
                        
end

%delete(gcp('nocreate'));
tt=toc;



%%
hold on; box on;

plot(FDT,squeeze(mean(average_sumSE(:,1,:),3)),'r-o','LineWidth',2);
plot(FDT,squeeze(mean(average_sumSE(:,2,:),3)),'b-<','LineWidth',2);
plot(FDT,squeeze(mean(average_sumSE2(:,:),2)),'m-+','LineWidth',2);

xlabel('Normalized Doppler Shift ($f_DT_s$)','Interpreter','latex');
ylabel('The sum SE (bit/s/Hz)','Interpreter','latex');
legend('aa=0.02','aa=0.05','by max FDT','Interpreter','latex');
set(gca, 'Fontname', 'Times New Roman','FontSize',14);
axis([0 0.004 0 60]);
grid on;




  

   

