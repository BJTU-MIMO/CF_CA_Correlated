


hold on; box on;
plot(FDT,squeeze(mean(average_sumSE2(:,:),2)),'r-o','LineWidth',2);
plot(FDT,squeeze(mean(average_sumSE1(:,:),2)),'b-+','LineWidth',2);
plot(FDT,squeeze(mean(average_sumSE0(:,:),2)),'k-<','LineWidth',2);


xlabel('Normalized Doppler Shift ($f_DT_s$)','Interpreter','latex');
ylabel('Sum SE (bit/s/Hz)','Interpreter','latex');
legend('By max $f_DT_s$','$\tau_c=200$','$\tau_c=1000$','Interpreter','latex');
set(gca, 'Fontname', 'Times New Roman','FontSize',14);
%axis([0 0.004 0 60]);
grid on;




  

   

