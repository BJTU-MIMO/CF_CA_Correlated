
clear;

parpool('local',4);
%% Define simulation setup

nbrOfSetups =200; %200;

nbrOfRealizations = 50;

L = 100;

K = 20;

N=2;

tau_c = 200;

tau_p =10;

tau_d=tau_c-tau_p;

data=0:1:(tau_d-1);

train=tau_p:-1:1;

Pmax = 100;

FDT=0.002;

rhotaup=besselj(0,2*pi*train*FDT);
rhotaup_hat=sqrt(1-rhotaup.^2);

rhon=besselj(0,2*pi*data*FDT);


SE_MR_simulation_tot_SC=zeros(nbrOfSetups,K);
SE_MR_simulation_tot1_SC=zeros(nbrOfSetups,K);

%%
parfor n = 1:nbrOfSetups    
     disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]); 
     
     [gainOverNoisedB,R,pilotIndex] = functionSetup(K,L,N,tau_p);
     gainOverNoise=db2pow(gainOverNoisedB);
             
     [beta,h,hhat,hc,gn,gamma_kl,Np] = functionChannelEstimates(R,gainOverNoise,K,L,N,nbrOfRealizations,Pmax,pilotIndex,rhotaup);     
         
     [SE_MR_simulation_SC] = functionComputeSE_AP_uplink_simulation_SC(hhat,R,gamma_kl,K,L,N,tau_c,tau_d,Pmax,rhon,nbrOfRealizations);
     SE_MR_simulation_tot_SC(n,:)=SE_MR_simulation_SC;
     SE_MR_simulation_tot1_SC(n)=mean(SE_MR_simulation_SC(:));



end
delete(gcp('nocreate')); 

%% Plot the simulation results
subplot(2,1,1);
hold on; box on;


plot(sort(reshape(SE_MR_simulation_tot_SC(:,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'r-','LineWidth',2);


xlabel('Per-User Uplink SE (bit/s/Hz)','Interpreter','latex');
ylabel('CDF','Interpreter','latex');
fig1=legend('small cell','Interpreter','latex');
set(gca, 'Fontname', 'Times New Roman','FontSize',14);
axis([0 5 0 1]);
grid on;

%save('smallcell0.mat','SE_MR_simulation_tot_SC');
%saveas(fig1,'smallcell0.fig');  


    


  

   

