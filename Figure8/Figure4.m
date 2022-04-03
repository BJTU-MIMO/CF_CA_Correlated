

hold on; box on;

plot(FDT,coherent_FPC(:,0.05*K*nbrOfSetups),'k-.<','LineWidth',2);
plot(FDT,coherent_taup(:,0.05*K*nbrOfSetups),'b--<','LineWidth',2);
plot(FDT,coherent(:,0.05*K*nbrOfSetups),'r-<','LineWidth',2);

plot(FDT,noncoherent_FPC(:,0.05*K*nbrOfSetups),'k-.o','LineWidth',2);
plot(FDT,noncoherent_taup(:,0.05*K*nbrOfSetups),'b--o','LineWidth',2);
plot(FDT,noncoherent(:,0.05*K*nbrOfSetups),'r-o','LineWidth',2);


xlabel('Normalized Doppler Shift ($f_DT_s$)','Interpreter','latex');
ylabel('95%-likely Per-User Downlink SE (bit/s/Hz)','Interpreter','latex');
legend('coherent (FPC)','coherent ($\tau_p=20$)','coherent ($\tau_p=10$)','non-coherent (FPC)','non-coherent ($\tau_p=20$)','non-coherent ($\tau_p=10$)','Interpreter','latex');
set(gca, 'Fontname', 'Times New Roman','FontSize',14);
axis([0 0.002 0 2]);
grid on;



  

   

