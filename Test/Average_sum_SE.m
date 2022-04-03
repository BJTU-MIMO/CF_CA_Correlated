
clear;

%% Define simulation setup

nbrOfSetups =10; %200;

nbrOfRealizations = 1000;

L = 20;

K = 10;

N=2;

tau_c = 200;

tau_p =5;

tau_d=tau_c-tau_p;

data=0:(tau_d-1);

train=tau_p:-1:1;

Pmax = 100;

FDT=0.002;

rhotaup=besselj(0,2*pi*train*FDT);
rhotaup_hat=sqrt(1-rhotaup.^2);

rhon=besselj(0,2*pi*FDT);

aa=zeros(nbrOfRealizations,K,L);
bb=zeros(K,L);

SE_MR_simulation_tot_CF=zeros(nbrOfSetups,K);
SE_MR_simulation_tot1_CF=zeros(nbrOfSetups);
SE_MR_analytical_tot_CF=zeros(nbrOfSetups,K);
SE_MR_analytical_tot1_CF=zeros(nbrOfSetups,K);

%%
for n = 1:nbrOfSetups    
     disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]); 
     
     [gainOverNoisedB,R,pilotIndex] = functionSetup(K,L,N,tau_p);
     gainOverNoise=db2pow(gainOverNoisedB);
             
     [beta,h,hhat,hc,gn,gamma_kl,Np] = functionChannelEstimates(R,gainOverNoise,K,L,N,nbrOfRealizations,Pmax,pilotIndex,rhotaup);
     
     [SE_MR_simulation_CF,term1,term2,term3,term4,term5] = functionComputeSE_AP_uplink_simulation_CF(hhat,h,gn,Np,beta,gamma_kl,K,L,tau_c,tau_d,Pmax,pilotIndex,rhon,nbrOfRealizations);
     SE_MR_simulation_tot_CF(n,:)=SE_MR_simulation_CF;
     SE_MR_simulation_tot1_CF(n)=mean(SE_MR_simulation_CF(:));
     
     [SE_MR_analytical_CF,terma1,terma2,terma3,terma4,terma5] = functionComputeSE_AP_uplink_analytical_CF(R,beta,gamma_kl,K,L,tau_c,tau_d,Pmax,pilotIndex,rhon);
     SE_MR_analytical_tot_CF(n,:)=SE_MR_analytical_CF;
     SE_MR_analytical_tot1_CF(n)=mean(SE_MR_analytical_CF(:));
    
%  for s=1:nbrOfRealizations
%          for k=1:K
%              for l=1:L
%                  aa(s,k,l)=hhat(:,s,k,l)'*h(:,s,k,l);
%                  bb(k,l)=trace(gamma_kl(:,:,k,l));
%              end
%          end
%      end
%      cc=squeeze(mean(aa(:,:,:),1));
% 
%      
           
end


%% Plot the simulation results

hold on; box on;

plot(sort(reshape(SE_MR_simulation_tot_CF(:,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'r-','LineWidth',2);
plot(sort(reshape(SE_MR_analytical_tot_CF(:,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'b-.','LineWidth',2);

xlabel('Per-User Uplink SE (bit/s/Hz)','Interpreter','latex');
ylabel('CDF','Interpreter','latex');
legend('Simulation','Analytical','Interpreter','latex');
set(gca, 'Fontname', 'Times New Roman','FontSize',14);
%axis([0 10000 0 1]);
grid on;



    


  

   

