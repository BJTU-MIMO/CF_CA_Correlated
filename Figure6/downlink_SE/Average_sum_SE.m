
clear;

%parpool('local',16);

%% Define simulation setup

nbrOfSetups =200; %200;

L = 100;

K = 20;

N=2;

tau_c = 200;

tau_p =10;

tau_d=tau_c-tau_p;

data=0:1:(tau_d-1);

train=tau_p:-1:1;

Pmax = 200;
P_pilot=100;

ASD=30;

FDT=[0 0.002];

SE_MR_analytical_tot_coherent=zeros(length(FDT),nbrOfSetups,K);
SE_MR_analytical_tot1_coherent=zeros(length(FDT),nbrOfSetups);

SE_MR_analytical_tot_noncoherent=zeros(length(FDT),nbrOfSetups,K);
SE_MR_analytical_tot1_noncoherent=zeros(length(FDT),nbrOfSetups);



%%
tic;
for n = 1:nbrOfSetups    
     disp(['Setup------------------- ' num2str(n) ' out of ' num2str(nbrOfSetups)]); 
     
         [gainOverNoisedB,R,pilotIndex] = functionSetup(K,L,N,tau_p,ASD);
               
         for s=1:length(FDT)
             disp(['FDT------ ' num2str(s) ' out of ' num2str(length(FDT))]); 
        
             rhotaup=besselj(0,2*pi*FDT(s)*train);
             rhon=besselj(0,2*pi*FDT(s)*data); 
     
             [gamma_kl,gamma_kil,tracQ] = functionChannelEstimates(R,K,L,N,P_pilot,pilotIndex,rhotaup);
         
            % full power 
            mu1=1./repmat(sum(tracQ(:,:),1),K,1);
                                
%%%%%%%%%%%     
     % full power
             [SE_MR_analytical_coherent] = functionComputeSE_AP_downlink_analytical_coherent(R,gamma_kl,gamma_kil,K,L,tau_c,tau_d,Pmax,mu1,pilotIndex,rhon);
             SE_MR_analytical_tot_coherent(s,n,:)=SE_MR_analytical_coherent;
             SE_MR_analytical_tot1_coherent(s,n)=mean(SE_MR_analytical_coherent(:));
 
             [SE_MR_analytical_noncoherent] = functionComputeSE_AP_downlink_analytical_noncoherent(R,gamma_kl,gamma_kil,K,L,tau_c,tau_d,Pmax,mu1,pilotIndex,rhon);
             SE_MR_analytical_tot_noncoherent(s,n,:)=SE_MR_analytical_noncoherent;
             SE_MR_analytical_tot1_noncoherent(s,n)=mean(SE_MR_analytical_noncoherent(:));
             
         end
               
end

%delete(gcp('nocreate'));
tt=toc;

%%
hold on; box on;

plot(sort(reshape(SE_MR_analytical_tot_coherent(1,:,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'r-','LineWidth',2);
plot(sort(reshape(SE_MR_analytical_tot_noncoherent(1,:,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'b--','LineWidth',2);
plot(sort(reshape(SE_MR_analytical_tot_coherent(2,:,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'m-.','LineWidth',2);
plot(sort(reshape(SE_MR_analytical_tot_noncoherent(2,:,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'k:','LineWidth',2.5);


xlabel('Per-User Downlink SE (bit/s/Hz)','Interpreter','latex');
ylabel('CDF','Interpreter','latex');
fig1=legend('Coherent (FDT=0)','Non-coherent (FDT=0)','Coherent (FDT=0.002)','Non-coherent (FDT=0.002)','Simulation','Interpreter','latex');
set(gca, 'Fontname', 'Times New Roman','FontSize',14);
%axis([0 4 0 1]);
grid on;


  

   

