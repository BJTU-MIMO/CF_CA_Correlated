function [SE] = functionComputeSE_AP_uplink_analytical_SCD(beta,gamma_kl,K,L,N,tau_c,tau_d,Pmax,pilotIndex,rho)

N0=1;
p=Pmax;

term1 = p*N^2*squeeze(abs(sum(gamma_kl(:,:),2)).^2);

test2 = zeros(L,K,K);
test3 = zeros(L,K,K);
for l=1:L
    for k=1:K
        for i=1:K
            test2(l,k,i) = gamma_kl(k,l)*beta(i,l);
            if i == k
                test3(l,k,i) = 0;
            else
                test3(l,k,i) = sqrt(gamma_kl(k,l)*gamma_kl(i,l));
            end

        end
    end
end
term2 = p*N*squeeze(sum(sum(test2(:,:,:),1),3));

term3=zeros(1,K);
for k=1:K
    term3(k)=p*N^2*sum(abs(sum(test3(:,k,(pilotIndex(k)==pilotIndex)'),1)).^2,3);
end

term4 = N0*N*squeeze(sum(gamma_kl(:,:),2));

SINR=zeros(K,tau_d);
for k=1:K
    for s=1:tau_d
        SINR(k,s)=rho(s)^2*term1(k)/(term2(k)+rho(s)^2*term3(k)+term4(k));
    end
end

SE=squeeze(sum(real(log2(1+SINR(:,:))),2))/tau_c;


