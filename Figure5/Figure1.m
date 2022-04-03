subplot(2,1,1);
hold on; box on;

plot(sort(reshape(SC0(:,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'k-.','LineWidth',2);
plot(sort(reshape(MF0(:,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'b--','LineWidth',2);
plot(sort(reshape(LSFD0(:,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'r-','LineWidth',2);

x1=[0.1 0.4 0.7 1.0 1.3 1.6 1.9 2.2 2.5 2.8 3.1 3.4 3.7 4.0 4.3 4.6 4.9 ...
    0.1 0.4 0.7 1.0 1.3 1.6 1.9 2.2 2.5 2.8 3.1 3.4 3.7 4.0];
y1=[0.0005 0.0034 0.011 0.027 0.060 0.12 0.21 0.34 0.49 0.64 0.77 0.87 0.94 0.97 0.99 1 1 ...
    0.0036 0.018 0.051 0.11 0.21 0.37 0.56 0.73 0.85 0.94 0.97 0.99 1 1];
plot(x1,y1,'ko','LineWidth',1);



xlabel('Per-User Uplink SE (bit/s/Hz)','Interpreter','latex');
ylabel('CDF','Interpreter','latex');
legend('SC','MF','LSFD','Simulation','Interpreter','latex');
set(gca, 'Fontname', 'Times New Roman','FontSize',14);
axis([0 5 0 1]);
grid on;

subplot(2,1,2);
hold on; box on;

plot(sort(reshape(SC02(:,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'k-.','LineWidth',2);
plot(sort(reshape(MF02(:,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'b--','LineWidth',2);
plot(sort(reshape(LSFD02(:,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'r-','LineWidth',2);

x2=[0.1 0.4 0.7 1.0 1.3 1.6 1.9 2.2 2.5 2.8 3.1 ...
    0.1 0.4 0.7 1.0 1.3 1.6 1.9 2.2 2.5];
y2=[0.001 0.015 0.059 0.18 0.36 0.58 0.77 0.90 0.96 0.99 1 ...
    0.007 0.059 0.22 0.49 0.75 0.91 0.98 1 1];
plot(x2,y2,'ko','LineWidth',1);

xlabel('Per-User Uplink SE (bit/s/Hz)','Interpreter','latex');
ylabel('CDF','Interpreter','latex');
legend('SC','MF','LSFD','Simulation','Interpreter','latex');
set(gca, 'Fontname', 'Times New Roman','FontSize',14);
axis([0 5 0 1]);
grid on;
