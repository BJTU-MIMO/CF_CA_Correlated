function [SE] = functionComputeSE_AP_uplink_analytical_SCD_limit(gamma_kl,K,L,tau_c,tau_d,Pmax,pilotIndex,rho)

p=Pmax;

term1 = p*squeeze(abs(sum(gamma_kl(:,:),2)).^2);

test3 = zeros(L,K,K);
for l=1:L
    for k=1:K
        for i=1:K
            if i == k
                test3(l,k,i) = 0;
            else
                test3(l,k,i) = sqrt(gamma_kl(k,l)*gamma_kl(i,l));
            end

        end
    end
end

term3=zeros(1,K);
for k=1:K
    term3(k)=p*sum(abs(sum(test3(:,k,(pilotIndex(k)==pilotIndex)'),1)).^2,3);
end

SINR=zeros(K,tau_d);
for k=1:K
    for s=1:tau_d
        SINR(k,s)=rho(s)^2*term1(k)/(rho(s)^2*term3(k));
    end
end

SE=squeeze(sum(real(log2(1+SINR(:,:))),2))/tau_c;


