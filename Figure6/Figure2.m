hold on; box on;

plot(sort(reshape(SE_MR_analytical_tot_coherent(1,:,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'r-','LineWidth',2);
plot(sort(reshape(SE_MR_analytical_tot_coherent(2,:,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'b--','LineWidth',2);
plot(sort(reshape(SE_MR_analytical_tot_noncoherent(1,:,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'m-.','LineWidth',2);
plot(sort(reshape(SE_MR_analytical_tot_noncoherent(2,:,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'k:','LineWidth',2.5);

x=[0.1 0.4 0.7 1.0 1.3 1.6 1.9 2.2 2.5 2.8 3.1 3.4 3.7 4.0 4.3 4.6 4.9 ...
    0.1 0.4 0.7 1.0 1.3 1.6 1.9 2.2 2.5 2.8 3.1 3.4 ...
    0.1 0.3 0.5 0.7 0.9 1.1 1.3 ...
    0.1 0.2 0.3 0.4 0.5 0.6 0.7];
y=[0.0026 0.013 0.030 0.057 0.11 0.20 0.33 0.48 0.64 0.78 0.89 0.95 0.98 0.99 1 1 1 ...
    0.006 0.032 0.11 0.27 0.50 0.72 0.88 0.95 0.99 1 1 1 ...
    0.014 0.089 0.24 0.44 0.67 0.88 0.99 ...
    0.044 0.16 0.32 0.53 0.74 0.91 0.99];
plot(x,y,'ko','LineWidth',1);

xlabel('Per-User Downlink SE (bit/s/Hz)','Interpreter','latex');
ylabel('CDF','Interpreter','latex');
fig1=legend('Coherent ($f_DT_s=0$)','Coherent ($f_DT_s=0.002$)','Non-coherent ($f_DT_s=0$)','Non-coherent ($f_DT_s=0.002$)','Simulation','Interpreter','latex');
set(gca, 'Fontname', 'Times New Roman','FontSize',14);
axis([0 5 0 1]);
grid on;