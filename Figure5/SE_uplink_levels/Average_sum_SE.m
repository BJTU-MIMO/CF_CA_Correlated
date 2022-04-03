
clear;

%parpool('local',16);

%% Define simulation setup

nbrOfSetups =200; %200;

L = 100;

K =20;

N=2;

tau_c = 200;

tau_p =10;

tau_d=tau_c-tau_p;

data=0:1:(tau_d-1);

train=tau_p:-1:1;

Pmax = 100;

ASD=30;

FDT=0.002;

SE_MR_analytical_tot_MF=zeros(nbrOfSetups,K);
SE_MR_analytical_tot1_MF=zeros(nbrOfSetups,1);

SE_MR_analytical_tot_LSFD=zeros(nbrOfSetups,K);
SE_MR_analytical_tot1_LSFD=zeros(nbrOfSetups,1);




%%
tic;
for n = 1:nbrOfSetups    
     disp(['Setup------------------- ' num2str(n) ' out of ' num2str(nbrOfSetups)]); 
     
         [gainOverNoisedB,R,pilotIndex] = functionSetup(K,L,N,tau_p,ASD);
        
     
             rhotaup=besselj(0,2*pi*FDT*train);
             rhon=besselj(0,2*pi*FDT*data); 
     
             [gamma_kl,gamma_kil] = functionChannelEstimates(R,K,L,N,Pmax,pilotIndex,rhotaup);
                          
%%%%%%%%%%%     
     % full power
             [SE_MR_analytical_MF] = functionComputeSE_AP_uplink_analytical_MF(R,gamma_kl,gamma_kil,K,L,tau_c,tau_d,Pmax,pilotIndex,rhon);
             SE_MR_analytical_tot_MF(n,:)=SE_MR_analytical_MF;
             SE_MR_analytical_tot1_MF(n)=mean(SE_MR_analytical_MF(:));
 
             [SE_MR_analytical_LSFD] = functionComputeSE_AP_uplink_analytical_LSFD(R,gamma_kl,gamma_kil,K,L,tau_c,tau_d,Pmax,pilotIndex,rhon);
             SE_MR_analytical_tot_LSFD(n,:)=SE_MR_analytical_LSFD;
             SE_MR_analytical_tot1_LSFD(n)=mean(SE_MR_analytical_LSFD(:));
             
             

end

%delete(gcp('nocreate'));
tt=toc;

%%

subplot(2,1,1);
hold on; box on;

plot(sort(reshape(SE_MR_analytical_tot_MF(:,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'b-','LineWidth',2);
plot(sort(reshape(SE_MR_analytical_tot_LSFD(:,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'r-','LineWidth',2);


xlabel('Per-User Uplink SE (bit/s/Hz)','Interpreter','latex');
ylabel('CDF','Interpreter','latex');
fig1=legend('MF','LSFD','Interpreter','latex');
set(gca, 'Fontname', 'Times New Roman','FontSize',14);
axis([0 5 0 1]);
grid on;


%% save

%save('CF02.mat','SE_MR_analytical_tot_MF','SE_MR_analytical_tot_LSFD');


  

   

