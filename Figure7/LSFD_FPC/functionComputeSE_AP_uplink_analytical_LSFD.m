function [SE] = functionComputeSE_AP_uplink_analytical_LSFD(R,gamma_kl,gamma_kil,K,L,tau_c,tau_d,Pmax,eta,pilotIndex,rho)

sigma2=1;
p=Pmax;

b_k=zeros(K,L);
Y=zeros(L,L,K,K);
c=zeros(L,K,K);
H_c=zeros(L,L,K,K);
A=zeros(L,L,K);

term2=zeros(L,L,K);
SINR=zeros(K,tau_d);

for l=1:L
    for k=1:K
        b_k(k,l)=trace(gamma_kl(:,:,k,l));
        A(l,l,k)=sigma2*trace(gamma_kl(:,:,k,l));
        for i=1:K
            Y(l,l,k,i)=p*eta(i)*trace(gamma_kl(:,:,k,l)*R(:,:,i,l));
            if i==k
               c(l,k,i)=0;
            else
               c(l,k,i)=trace(gamma_kil(:,:,k,i,l)); 
            end            
        end
    end
end
term1=squeeze(sum(Y(:,:,:,:),4));

for k=1:K
    for i=1:K
    H_c(:,:,k,i)=p*eta(i)*c(:,k,i)*c(:,k,i)';    
    end
end

for k=1:K
    term2(:,:,k)=sum(H_c(:,:,k,(pilotIndex(k)==pilotIndex)'),4);
end
term3=A(:,:,:);

for k=1:K
    for s=1:tau_d
        SINR(k,s)=rho(s)^2*p*eta(k)*b_k(k,:)*(term1(:,:,k)+rho(s)^2*term2(:,:,k)+term3(:,:,k))^(-1)*b_k(k,:)';
    end
end

SE=squeeze(sum(real(log2(1+SINR(:,:))),2))/tau_c;



 
