
clear;

%parpool('local',16);

%% Define simulation setup

nbrOfSetups =200; %200;

L = 100;

K = 20;

N=2;

tau_c = 500;

tau_p =10;

tau_d=tau_c-tau_p;

data_u=0:1:(tau_d-1);
data_d=0:1:(tau_d-1);

train=tau_p:-1:1;

P_pilot = 100;
P_u=100;
P_d=200;

ASD=30;

FDT=[0.001 0.002];

SE_MR_analytical_tot_LSFD=zeros(length(FDT),nbrOfSetups,K,tau_d);
SE_MR_analytical_tot1_LSFD=zeros(length(FDT),nbrOfSetups,tau_d);

SE_MR_analytical_tot_coherent=zeros(length(FDT),nbrOfSetups,K,tau_d);
SE_MR_analytical_tot1_coherent=zeros(length(FDT),nbrOfSetups,tau_d);


%%
tic;
for n = 1:nbrOfSetups    
     disp(['Setup------------------- ' num2str(n) ' out of ' num2str(nbrOfSetups)]); 
     
         [gainOverNoisedB,R,pilotIndex] = functionSetup(K,L,N,tau_p,ASD);
               
         for s=1:length(FDT)
             disp(['FDT------ ' num2str(s) ' out of ' num2str(length(FDT))]); 
        
             rhotaup=besselj(0,2*pi*FDT(s)*train);
             rhon_u=besselj(0,2*pi*FDT(s)*data_u); 
             rhon_d=besselj(0,2*pi*FDT(s)*data_d); 
     
             [gamma_kl,gamma_kil,tracQ] = functionChannelEstimates(R,K,L,N,P_pilot,pilotIndex,rhotaup);
         
            % full power 
              mu1=1./repmat(sum(tracQ(:,:),1),K,1);
                          
%%%%%%%%%%%     

             [SE_MR_analytical_LSFD] = functionComputeSE_AP_uplink_analytical_LSFD(R,gamma_kl,gamma_kil,K,L,tau_d,P_u,pilotIndex,rhon_u);
             SE_MR_analytical_tot_LSFD(s,n,:,:)=SE_MR_analytical_LSFD;
             SE_MR_analytical_tot1_LSFD(s,n,:)=squeeze(mean(SE_MR_analytical_LSFD(:,:),1));
     
             [SE_MR_analytical_coherent] = functionComputeSE_AP_downlink_analytical_coherent(R,gamma_kl,gamma_kil,K,L,tau_d,P_d,mu1,pilotIndex,rhon_d);
             SE_MR_analytical_tot_coherent(s,n,:,:)=SE_MR_analytical_coherent;
             SE_MR_analytical_tot1_coherent(s,n,:)=squeeze(mean(SE_MR_analytical_coherent(:,:),1));
             
        end
    
end

%delete(gcp('nocreate'));
tt=toc;

%%
hold on; box on;

plot(data_u+11,squeeze(mean(SE_MR_analytical_tot1_LSFD(1,:,:),2)),'r-','LineWidth',2);
plot(data_d+11,squeeze(mean(SE_MR_analytical_tot1_coherent(1,:,:),2)),'b--','LineWidth',2);
plot(data_u+11,squeeze(mean(SE_MR_analytical_tot1_LSFD(2,:,:),2)),'m-.','LineWidth',2);
plot(data_d+11,squeeze(mean(SE_MR_analytical_tot1_coherent(2,:,:),2)),'k:','LineWidth',2.5);

xlabel('Time Instant Index $n$','Interpreter','latex');
ylabel('Average SE$[n]$','Interpreter','latex');
fig1=legend('LSFD (FDT=0.001)','coherent (FDT=0.001)','LSFD (FDT=0.002)','coherent (FDT=0.002)','Interpreter','latex');
set(gca, 'Fontname', 'Times New Roman','FontSize',14);
%axis([0 500 0 2]);
grid on;
annotation('textarrow',[0.24 0.15],[0.88 0.84],'String','$n=\tau_p+1$','Interpreter','latex');
annotation('textarrow',[0.20 0.15],[0.58 0.72],'String','$n=\tau_p+2$','Interpreter','latex');

%% save

% save('CF.mat','tt','SE_MR_analytical_tot_MF','SE_MR_analytical_tot_LSFD','SE_MR_analytical_tot_FPC_MF','SE_MR_analytical_tot_FPC_LSFD','SE_MR_analytical_tot_maxmin_MF');
% saveas(fig1,'CF.fig');  

  

   

