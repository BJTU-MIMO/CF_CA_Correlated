function [SE] = functionComputeSE_AP_downlink_analytical_noncoherent(R,gamma_kl,gamma_kil,K,L,tau_c,tau_d,Pmax,mu,pilotIndex,rho)

sigma2=1;
p=Pmax;

test1=zeros(K,L);
test2=zeros(K,K,L);
test3=zeros(K,K,L);
term3=zeros(1,K);
SINR=zeros(K,tau_d);

for l=1:L
    for k=1:K
        test1(k,l)=sqrt(mu(k,l))*trace(gamma_kl(:,:,k,l));
        for i=1:K
            test2(k,i,l)=mu(i,l)*trace(gamma_kl(:,:,i,l)*R(:,:,k,l));
            if i==k
                test3(k,i,l)=0;
            else
                test3(k,i,l)=mu(i,l)*trace(gamma_kil(:,:,k,i,l))*trace(gamma_kil(:,:,k,i,l));
            end
        end
    end
end
term1=squeeze(p*sum(abs(test1(:,:)).^2,2));
term2=squeeze(p*sum(sum(test2(:,:,:),3),2));
test33=squeeze(sum(test3(:,:,:),3));

for k=1:K
    term3(k)=p*sum(test33(k,(pilotIndex(k)==pilotIndex)'),2);
end

for k=1:K
    for s=1:tau_d
        SINR(k,s)=rho(s)^2*term1(k)/(term2(k)+rho(s)^2*term3(k)+sigma2);
    end
end

SE=squeeze(sum(real(log2(1+SINR(:,:))),2))/tau_c;
 
