
hold on; box on;

plot(L,squeeze(mean(EE_tot_K40(:,1,:),3)),'r-o','LineWidth',2);
plot(L,squeeze(mean(EE_tot_K20(:,1,:),3)),'b--+','LineWidth',2);

plot(L,squeeze(mean(EE_tot_K40(:,2,:),3)),'m-.x','LineWidth',2);
plot(L,squeeze(mean(EE_tot_K20(:,2,:),3)),'k:<','LineWidth',2.5);


xlabel('Number of APs ($L$)','Interpreter','latex');
ylabel('Total EE (Mbit/s/Hz)','Interpreter','latex');
legend('$f_DT_s=0.001$, $K=40$','$f_DT_s=0.001$, $K=20$','$f_DT_s=0.002$, $K=40$','$f_DT_s=0.002$, $K=20$','Interpreter','latex');
set(gca, 'Fontname', 'Times New Roman','FontSize',14);
axis([10 100 4 12]);
grid on;

%%


hold on; box on;

plot(L,EE_K40(:,1),'r-o','LineWidth',2);
plot(L,EE_K20(:,1),'b--+','LineWidth',2);

plot(L,EE_K40(:,2),'m-.x','LineWidth',2);
plot(L,EE_K20(:,2),'k:<','LineWidth',2.5);


xlabel('Number of APs ($L$)','Interpreter','latex');
ylabel('Total EE (Mbit/s/Hz)','Interpreter','latex');
legend('$f_DT_s=0.001$, $K=40$','$f_DT_s=0.001$, $K=20$','$f_DT_s=0.002$, $K=40$','$f_DT_s=0.002$, $K=20$','Interpreter','latex');
set(gca, 'Fontname', 'Times New Roman','FontSize',14);
axis([10 100 4 12]);
grid on;



