
hold on; box on;

plot(sort(reshape(SE_MR_analytical_tot_coherent(1,1,:,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'r-','LineWidth',2);
plot(sort(reshape(SE_MR_analytical_tot_coherent(1,2,:,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'b-','LineWidth',2);
plot(sort(reshape(SE_MR_analytical_tot_coherent(2,1,:,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'r--','LineWidth',2);
plot(sort(reshape(SE_MR_analytical_tot_coherent(2,2,:,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'b--','LineWidth',2);
plot(sort(reshape(SE_MR_analytical_tot_coherent_UN(1,1,:,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'r-.','LineWidth',2);
plot(sort(reshape(SE_MR_analytical_tot_coherent_UN(1,2,:,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'b-.','LineWidth',2);




xlabel('Per-User Downlink SE (bit/s/Hz)','Interpreter','latex');
ylabel('CDF','Interpreter','latex');
legend('$\mathrm{ASD}=10^\mathrm{o}$, $f_DT_s=0.001$','$\mathrm{ASD}=10^\mathrm{o}$, $f_DT_s=0.002$','$\mathrm{ASD}=50^\mathrm{o}$, $f_DT_s=0.001$','$\mathrm{ASD}=50^\mathrm{o}$, $f_DT_s=0.002$','uncorrelated, $f_DT_s=0.001$','uncorrelated, $f_DT_s=0.002$','Interpreter','latex');
set(gca, 'Fontname', 'Times New Roman','FontSize',14);
%axis([10 50 0 4]);
grid on;


%% save

%save('CF02.mat','SE_MR_analytical_tot_MF','SE_MR_analytical_tot_LSFD');


  

   

