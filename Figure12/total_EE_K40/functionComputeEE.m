function [EE] = functionComputeEE(SE_u,SE_d,P_u,P_d,tracQ,eta,mu1,partial,P_ue,P_ap,P_0,P_bt,N,K,L,B,tau_c,tau_p)

sigma2=1;
P_u=P_u/1000;
P_d=P_d/1000;

SE_sum=sum(0.5*SE_u(:)+0.5*SE_d(:));
P1=K*(((tau_c+tau_p)/(2*tau_c))*P_u*sigma2*eta/partial+P_ue);
P2=sum(((tau_c-tau_p)/(2*tau_c))*P_d*sigma2/partial*squeeze(sum(mu1(:,:).*tracQ(:,:),1))+N*P_ap+P_0);
P3=B*L*P_bt*SE_sum/(10e9);
P_total=P1+P2+P3;
EE=B*SE_sum/P_total;