
clear;

%parpool('local',4);
%% Define simulation setup

nbrOfSetups =200; %200;

nbrOfRealizations = 50;

L = 100;

K = 20;

N=2;

tau_c = 200;

tau_p =20;

tau_d=tau_c-tau_p;

data=0:1:(tau_d-1);

train=tau_p:-1:1;

Pmax = 100;

FDT = 0:0.0005:0.003;


SE_MR_simulation_tot_SC=zeros(length(FDT),nbrOfSetups,K);
SE_MR_simulation_tot1_SC=zeros(length(FDT),nbrOfSetups);

%%
for n = 1:nbrOfSetups    
     disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]); 
     
     [gainOverNoisedB,R,pilotIndex] = functionSetup(K,L,N,tau_p);
     gainOverNoise=db2pow(gainOverNoisedB);
     
     for s=1:length(FDT)
         disp(['FDT ' num2str(s) ' out of ' num2str(length(FDT))]); 
     
     rhotaup=besselj(0,2*pi*train*FDT(s));
     rhon=besselj(0,2*pi*data*FDT(s));
             
     [beta,h,hhat,hc,gn,gamma_kl,Np] = functionChannelEstimates(R,gainOverNoise,K,L,N,nbrOfRealizations,Pmax,pilotIndex,rhotaup);     
         
     [SE_MR_simulation_SC] = functionComputeSE_AP_uplink_simulation_SC(hhat,R,gamma_kl,K,L,N,tau_c,tau_d,Pmax,rhon,nbrOfRealizations);
     SE_MR_simulation_tot_SC(s,n,:)=SE_MR_simulation_SC;
     SE_MR_simulation_tot1_SC(s,n)=mean(SE_MR_simulation_SC(:));
     
     end

end
%delete(gcp('nocreate')); 


%%
smallcell_taup=zeros(length(FDT),K*nbrOfSetups);

for i=1:length(FDT)
smallcell_taup(i,:)=sort(reshape(SE_MR_simulation_tot_SC(i,:,:),[K*nbrOfSetups 1]));
end


%%

hold on; box on;

plot(FDT,smallcell_taup(:,0.05*K*nbrOfSetups),'r-o','LineWidth',2);

xlabel('Normalized Doppler Shift ($f_DT_s$)','Interpreter','latex');
ylabel('95%-likely Per-User Uplink SE (bit/s/Hz)','Interpreter','latex');
legend('smallcell_taup','Interpreter','latex');
set(gca, 'Fontname', 'Times New Roman','FontSize',14);
%axis([0 4 0 1]);
grid on;


    


  

   

