function [SE,term1,term2,term3,term4,term5] = functionComputeSE_AP_uplink_analytical_CF(R,beta,gamma_kl,K,L,tau_c,tau_d,Pmax,pilotIndex,rho)

sigma2=1;
p=Pmax;

test1=zeros(K,L);
test2=zeros(K,L);
test3=zeros(K,K,L);
test4=zeros(K,K,L);

term4=zeros(1,K);

for k=1:K
    for l=1:L
        test1(k,l)=trace(gamma_kl(:,:,k,l));
        test2(k,l)=trace(gamma_kl(:,:,k,l)*R(:,:,k,l));
        for i=1:K
            if i==k
            test3(k,i,l)=0;
            test4(k,i,l)=0;
            else
            test3(k,i,l)=trace(gamma_kl(:,:,k,l)*R(:,:,i,l));
            test4(k,i,l)=sqrt(trace(gamma_kl(:,:,k,l))*trace(gamma_kl(:,:,i,l)));
            end
        end
    end
end
term1=rho^2*p*abs(squeeze(sum(test1(:,:),2))).^2;
term2=rho^2*p*squeeze(sum(test2(:,:),2));
term3=(1-rho^2)*p*squeeze(sum(test2(:,:),2));
for k=1:K
    term4(k)=sum(p*sum(test3(k,:,:),3),2)+sum(rho^2*p*abs(sum(test4(k,(pilotIndex(k)==pilotIndex)',:),3)).^2,2);
end
term5=sigma2*squeeze(sum(abs(test1(:,:)),2));

SINR=term1(:)./(term2(:)+term3(:)+term4(:)+term5(:));
SE=real(log2(1+SINR(:)));




