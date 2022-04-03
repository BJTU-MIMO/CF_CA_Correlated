
clear;

%parpool('local',16);

%% Define simulation setup

nbrOfSetups =5; %200;

Bandwith=20e6;

%L = [10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100];% 50 60 70 80 90 100];
L = 10:10:100;

K =20;

N=2;

tau_c = 200;

tau_p =10;

tau_d=tau_c-tau_p;

data_u=0:1:(tau_d-1);
data_d=0:1:(tau_d-1);

train=tau_p:-1:1;

P_pilot = 100;
P_u=100;
P_d=200;

eta=1;

partial=0.4;
P_ue=0.1;
P_ap=0.2;
P_0=0.825;
P_bt=0.25;

ASD=10;

FDT=[0.001 0.002];

SE_MR_analytical_tot_LSFD=zeros(length(L),length(FDT),nbrOfSetups,K);
SE_MR_analytical_tot1_LSFD=zeros(length(L),length(FDT),nbrOfSetups);
SE_MR_analytical_tot_coherent=zeros(length(L),length(FDT),nbrOfSetups,K);
SE_MR_analytical_tot1_coherent=zeros(length(L),length(FDT),nbrOfSetups);

EE_tot=zeros(length(L),length(FDT),nbrOfSetups);

%%
tic;
for n = 1:nbrOfSetups    
     disp(['Setup------------------- ' num2str(n) ' out of ' num2str(nbrOfSetups)]); 
     for r=1:length(L)
         disp(['L----' num2str(r) ' out of ' num2str(length(L))]); 
         for s=1:length(FDT)
             disp(['FDT------ ' num2str(s) ' out of ' num2str(length(FDT))]); 
     
             [gainOverNoisedB,R,pilotIndex] = functionSetup(K,L(r),N,tau_p,ASD,Bandwith);

                     
             rhotaup=besselj(0,2*pi*FDT(s)*train);
             rhon_u=besselj(0,2*pi*FDT(s)*data_u); 
             rhon_d=besselj(0,2*pi*FDT(s)*data_d); 
     
             [gamma_kl,gamma_kil,tracQ] = functionChannelEstimates(R,K,L(r),N,P_pilot,pilotIndex,rhotaup);
         
            % full power
             mu1=1./repmat(sum(tracQ(:,:),1),K,1);  
                          
%%%%%%%%%%%     
     % full power
 
             [SE_MR_analytical_LSFD] = functionComputeSE_AP_uplink_analytical_LSFD(R,gamma_kl,gamma_kil,K,L(r),tau_c,tau_d,P_u,pilotIndex,rhon_u);
             SE_MR_analytical_tot_LSFD(r,s,n,:)=SE_MR_analytical_LSFD;
             SE_MR_analytical_tot1_LSFD(r,s,n)=mean(SE_MR_analytical_LSFD(:));
             
             [SE_MR_analytical_coherent] = functionComputeSE_AP_downlink_analytical_coherent(R,gamma_kl,gamma_kil,K,L(r),tau_c,tau_d,P_d,mu1,pilotIndex,rhon_d);
             SE_MR_analytical_tot_coherent(r,s,n,:)=SE_MR_analytical_coherent;
             SE_MR_analytical_tot1_coherent(r,s,n)=mean(SE_MR_analytical_coherent(:));
             
             [EE] = functionComputeEE(SE_MR_analytical_LSFD,SE_MR_analytical_coherent,P_u,P_d,tracQ,eta,mu1,partial,P_ue,P_ap,P_0,P_bt,N,K,L(r),Bandwith,tau_c,tau_p);
             EE_tot(r,s,n)=EE;
             
         end
     end
                         
end

%delete(gcp('nocreate'));
tt=toc;

%%
%LSFD=squeeze(mean(SE_MR_analytical_tot1_LSFD(:,:,:),3));



%%
EE_tot=EE_tot./1000000;

%%
hold on; box on;

plot(L,squeeze(mean(EE_tot(:,1,:),3)),'r-','LineWidth',2);
plot(L,squeeze(mean(EE_tot(:,2,:),3)),'b-','LineWidth',2);

xlabel('The number of APs ($L$)','Interpreter','latex');
ylabel('The total energy efficiency (Mbit/s/Hz)','Interpreter','latex');
legend('FDT=0','FDT=0.003','Interpreter','latex');
set(gca, 'Fontname', 'Times New Roman','FontSize',14);
%axis([10 100 0 10]);
grid on;


%% save

%save('CF02.mat','SE_MR_analytical_tot_MF','SE_MR_analytical_tot_LSFD');


  

   

