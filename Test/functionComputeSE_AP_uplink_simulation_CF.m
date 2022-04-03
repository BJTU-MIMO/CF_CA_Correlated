function [SE,term1,term2,term3,term4,term5] = functionComputeSE_AP_uplink_simulation_CF(hhat,h,gn,Np,beta,gamma_kl,K,L,tau_c,tau_d,Pmax,pilotIndex,rho,nbrOfRealizations)

sigma2=1;
p=Pmax;

test1=zeros(nbrOfRealizations,K,L);
test2=zeros(nbrOfRealizations,K,L);
test3=zeros(nbrOfRealizations,K,L);
test4=zeros(nbrOfRealizations,K,K,L);

     for s=1:nbrOfRealizations
         for k=1:K
             for l=1:L
                 test1(s,k,l)=hhat(:,s,k,l)'*h(:,s,k,l);
                 test2(s,k,l)=hhat(:,s,k,l)'*gn(:,s,k,l);
                 test3(s,k,l)=hhat(:,s,k,l)'*Np(:,s,k,l);
                 for i=1:K
                     if i==k
                        test4(s,k,i,l)=0;
                     else
                        test4(s,k,i,l)=hhat(:,s,k,l)'*h(:,s,i,l);
                     end
                 end
                 
             end
         end
     end

term1=rho^2*p*abs(squeeze(sum(squeeze(mean(test1(:,:,:),1)),2))).^2;
term2=rho^2*p*squeeze(sum((mean(abs(test1(:,:,:)).^2,1)-abs(mean(test1(:,:,:),1)).^2),3));
term3=(1-rho^2)*p*squeeze(mean(abs(sum(test2(:,:,:),3)).^2,1));
term4=squeeze(sum(p*mean(abs(sum(test4(:,:,:,:),4)).^2,1),3));
term5=squeeze(mean(abs(sum(test3(:,:,:),3)).^2,1));

SINR=term1(:)./(term2(:)+term3(:)+term4(:)+term5(:));
SE=real(log2(1+SINR(:)));




