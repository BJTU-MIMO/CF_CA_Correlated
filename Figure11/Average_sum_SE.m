
clear;

%% Define simulation setup

nbrOfSetups =200; %200;

nbrOfRealizations = 1000;

L = 100;

K = 20;

tau_c = 200;

tau_p =[10 20];

Pmax = 100;

FDT=0.002;

%N = 10:10:100;
N =1:15;
%N = [1:9 10:10:90 100:100:900 1000:1000:50000];
%N = [1000 2000 3000 4000 5000 7000 8000 9000 10000 20000 30000 40000 50000];
%N = [10 20 30 40 50 60 70 80 90 100 150 200 250 300 350 400 450 500];

SE_MR_analytical_tot_SCD=zeros(length(N),length(tau_p),nbrOfSetups,K);
SE_MR_analytical_tot1_SCD=zeros(length(N),length(tau_p),nbrOfSetups);

SE_MR_analytical_tot_SCD_limit=zeros(length(N),nbrOfSetups,K);
SE_MR_analytical_tot1_SCD_limit=zeros(length(N),nbrOfSetups);


%%
for s = 1:nbrOfSetups
    disp(['Setup ' num2str(s) ' out of ' num2str(nbrOfSetups)]);
    
    for t = 1:length(tau_p)
        
        [gainOverNoisedB,pilotIndex] = functionSetup(K,L,tau_p(t));
        gainOverNoise=db2pow(gainOverNoisedB);
        
        tau_d=tau_c-tau_p(t);
        data=0:(tau_d-1);
        train=tau_p(t):-1:1;
        
        rhon=besselj(0,2*pi*FDT*data);
        rhotaup=besselj(0,2*pi*train*FDT);
        rhotaup_hat=sqrt(1-rhotaup.^2);
        
        [beta,h,hhat,hc,gn,gamma_kl,Np] = functionChannelEstimates(gainOverNoise,K,L,nbrOfRealizations,Pmax,pilotIndex,rhotaup,rhotaup_hat);
        
        for n = 1:length(N)
            
            [SE_MR_analytical_SCD] = functionComputeSE_AP_uplink_analytical_SCD(beta,gamma_kl,K,L,N(n),tau_c,tau_d,Pmax,pilotIndex,rhon);
            SE_MR_analytical_tot_SCD(n,t,s,:)=SE_MR_analytical_SCD;
            SE_MR_analytical_tot1_SCD(n,t,s)=mean(SE_MR_analytical_SCD(:));
            
            if t == 1
                [SE_MR_analytical_SCD_limit] = functionComputeSE_AP_uplink_analytical_SCD_limit(gamma_kl,K,L,tau_c,tau_d,Pmax,pilotIndex,rhon);
                SE_MR_analytical_tot_SCD_limit(n,s,:)=SE_MR_analytical_SCD_limit;
                SE_MR_analytical_tot1_SCD_limit(n,s)=mean(SE_MR_analytical_SCD_limit(:));
            end
            
        end
        
    end
    
end



%% Plot the simulation results

% hold on; box on;
% 
% plot(N,squeeze(mean(SE_MR_analytical_tot1_SCD(:,1,:),3)),'r-','LineWidth',2);
% plot(N,squeeze(mean(SE_MR_analytical_tot1_SCD(:,2,:),3)),'b-','LineWidth',2);
% plot(N,squeeze(mean(SE_MR_analytical_tot1_SCD_limit(:,:),2)),'k--','LineWidth',2);
% 
% xlabel('Number of antennas per APs ($N$)','Interpreter','latex');
% ylabel('Average Uplink SE (bit/s/Hz)','Interpreter','latex');
% legend('$\tau_p=10$ (pilot contamination)','$\tau_p=20$ (no pilot contamination)','Limit with $\tau_p=10$','Interpreter','latex');
% set(gca, 'Fontname', 'Times New Roman','FontSize',14);
% axis([1 50000 0 14]);
% grid on;

%%
%magnify;

%%

hold on; box on;

plot(N,squeeze(mean(SE_MR_analytical_tot1_SCD(:,2,:),3)),'r-o','LineWidth',2);
plot(N,squeeze(mean(SE_MR_analytical_tot1_SCD(:,1,:),3)),'b-->','LineWidth',2);
%plot(N,squeeze(mean(SE_MR_analytical_tot1_SCD_limit(:,:),2)),'k--','LineWidth',2);

xlabel('Number of antennas per APs ($N$)','Interpreter','latex');
ylabel('Average Uplink SE (bit/s/Hz)','Interpreter','latex');
legend('$\tau_p=20$ (no pilot contamination)','$\tau_p=10$ (pilot contamination)','Interpreter','latex');
set(gca, 'Fontname', 'Times New Roman','FontSize',14);
axis([1 15 0 3.5]);
grid on;








