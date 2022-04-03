function [SE] = functionComputeSE_AP_uplink_simulation_SC(hhat,R,gamma_kl,K,L,N,tau_c,tau_d,Pmax,rho,nbrOfRealizations)

sigma2=1;
p=Pmax;

test1=zeros(nbrOfRealizations,K,L);
test2=zeros(nbrOfRealizations,K,K,L);
test3=zeros(nbrOfRealizations,K,K,L);
test4=zeros(nbrOfRealizations,K,K,L);
term5=zeros(nbrOfRealizations,K,L);

SINR=zeros(nbrOfRealizations,K,L,tau_d);

     for s=1:nbrOfRealizations
         for k=1:K
             for l=1:L
                 test1(s,k,l)=abs(hhat(:,s,k,l)'*hhat(:,s,k,l))^2;
                 term5(s,k,l)=sigma2*hhat(:,s,k,l)'*eye(N)*hhat(:,s,k,l);
                 for i=1:K
                     test3(s,k,i,l)=hhat(:,s,k,l)'*R(:,:,i,l)*hhat(:,s,k,l);
                     test4(s,k,i,l)=hhat(:,s,k,l)'*gamma_kl(:,:,i,l)*hhat(:,s,k,l);
                     if i==k
                        test2(s,k,i,l)=0;
                     else
                        test2(s,k,i,l)=abs(hhat(:,s,k,l)'*hhat(:,s,i,l))^2;
                     end
                 end
                 
             end
         end
     end

term1=p*test1(:,:,:);
term2=p*squeeze(sum(test2(:,:,:,:),3));
term3=p*squeeze(sum(test3(:,:,:,:),3));
term4=p*squeeze(sum(test4(:,:,:,:),3));

for t=1:tau_d   
    SINR(:,:,:,t)=rho(t)^2*term1(:,:,:)./(rho(t)^2*term2(:,:,:)+term3(:,:,:)-rho(t)^2*term4(:,:,:)+term5(:,:,:));
end

SE_kl=squeeze(sum(squeeze(mean(real(log2(1+SINR(:,:,:,:))),1)),3))/tau_c;
SE=max(SE_kl,[],2);





