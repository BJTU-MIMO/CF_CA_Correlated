
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

FDT=0:0.0005:0.003;


SE_MR_analytical_tot_LSFD=zeros(length(FDT),nbrOfSetups,K);
SE_MR_analytical_tot1_LSFD=zeros(length(FDT),nbrOfSetups);


%%
tic;
for n = 1:nbrOfSetups    
     disp(['Setup------------------- ' num2str(n) ' out of ' num2str(nbrOfSetups)]); 
     
         [gainOverNoisedB,R,pilotIndex] = functionSetup(K,L,N,tau_p,ASD);
         
         %fractional power control
         beta=db2pow(gainOverNoisedB);
         beta_sum=squeeze(sum(beta(:,:),2));
         beta_min=min(beta_sum);
         beta_max=max(beta_sum);        
         eta_FPC=beta_min./beta_sum;  %差用户分的功率多
          %eta_FPC=ones(1,K);
         %eta_FPC=beta_sum/beta_max;   %好用户分的多
         
         for s=1:length(FDT)
             disp(['FDT------ ' num2str(s) ' out of ' num2str(length(FDT))]); 
             
             rhotaup=besselj(0,2*pi*FDT(s)*train);
             rhon=besselj(0,2*pi*FDT(s)*data); 
     
             [gamma_kl,gamma_kil] = functionChannelEstimates(R,K,L,N,Pmax,pilotIndex,rhotaup);
                          
%%%%%%%%%%%     
     % full power
 
             [SE_MR_analytical_LSFD] = functionComputeSE_AP_uplink_analytical_LSFD(R,gamma_kl,gamma_kil,K,L,tau_c,tau_d,Pmax,eta_FPC,pilotIndex,rhon);
             SE_MR_analytical_tot_LSFD(s,n,:)=SE_MR_analytical_LSFD;
             SE_MR_analytical_tot1_LSFD(s,n)=mean(SE_MR_analytical_LSFD(:));
             
         end
                         
end

%delete(gcp('nocreate'));
tt=toc;

%%
LSFD_FPC=zeros(length(FDT),K*nbrOfSetups);

for i=1:length(FDT)
LSFD_FPC(i,:)=sort(reshape(SE_MR_analytical_tot_LSFD(i,:,:),[K*nbrOfSetups 1]));

end


%%

hold on; box on;

plot(FDT,LSFD_FPC(:,0.05*K*nbrOfSetups),'b-o','LineWidth',2);

xlabel('Normalized Doppler Shift ($f_DT_s$)','Interpreter','latex');
ylabel('95%-likely Per-User Uplink SE (bit/s/Hz)','Interpreter','latex');
legend('LSFD_FPC','Interpreter','latex');
set(gca, 'Fontname', 'Times New Roman','FontSize',14);
%axis([0 4 0 1]);
grid on;


%% save

%save('CF02.mat','SE_MR_analytical_tot_MF','SE_MR_analytical_tot_LSFD');


  

   

