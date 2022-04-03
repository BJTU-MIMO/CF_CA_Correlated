

hold on; box on;
plot(FDT,LSFD_FPC(:,0.05*K*nbrOfSetups),'k-.o','LineWidth',2);
plot(FDT,LSFD_taup(:,0.05*K*nbrOfSetups),'b--o','LineWidth',2);
plot(FDT,LSFD(:,0.05*K*nbrOfSetups),'r-o','LineWidth',2);

plot(FDT,smallcell_FPC(:,0.05*K*nbrOfSetups),'k-.<','LineWidth',2);
plot(FDT,smallcell_taup(:,0.05*K*nbrOfSetups),'b--<','LineWidth',2);
plot(FDT,smallcell(:,0.05*K*nbrOfSetups),'r-<','LineWidth',2);



xlabel('Normalized Doppler Shift ($f_DT_s$)','Interpreter','latex');
ylabel('95%-likely Per-User Uplink SE (bit/s/Hz)','Interpreter','latex');
legend('LSFD (FPC, $\tau_{p}=10$)','LSFD (full power, $\tau_{p}=20$)','LSFD (full power, $\tau_{p}=10$)','SC (FPC, $\tau_{p}=10$)','SC (full power, $\tau_{p}=20$)','SC (full power, $\tau_{p}=10$)','Interpreter','latex');
set(gca, 'Fontname', 'Times New Roman','FontSize',14);
axis([0 0.002 0 2.2]);
grid on;




  

   

