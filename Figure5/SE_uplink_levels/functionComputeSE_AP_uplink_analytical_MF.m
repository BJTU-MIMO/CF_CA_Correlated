function [SE] = functionComputeSE_AP_uplink_analytical_MF(R,gamma_kl,gamma_kil,K,L,tau_c,tau_d,Pmax,pilotIndex,rho)

sigma2=1;
p=Pmax;

a_kl=ones(K,L);

b_k=zeros(K,L);
Y=zeros(L,L,K,K);
c_ki=zeros(K,K,L);
A=zeros(L,L,K);

term1=zeros(1,K);
term3=zeros(1,K);
term4=zeros(1,K);

test2=zeros(K,K);
test3=zeros(K,K);

SINR=zeros(K,tau_d);

a_k=a_kl(:,:);

for l=1:L
    for k=1:K
        b_k(k,l)=trace(gamma_kl(:,:,k,l));
        A(l,l,k)=trace(gamma_kl(:,:,k,l));
        for i=1:K
            Y(l,l,k,i)=trace(gamma_kl(:,:,k,l)*R(:,:,i,l));
            c_ki(k,i,l)=trace(gamma_kil(:,:,k,i,l));
        end
    end
end

for k=1:K
    term1(k)=p*abs(a_k(k,:)*b_k(k,:)')^2;
    term4(k)=sigma2*a_k(k,:)*A(:,:,k)*a_k(k,:)';
    for i=1:K
        test2(k,i)=p*a_k(k,:)*Y(:,:,k,i)*a_k(k,:)';
        if i==k
           test3(k,i)=0;
        else
           test3(k,i)=p*abs(a_k(k,:)*squeeze(c_ki(k,i,:)))^2; 
        end
    end
end

term2=squeeze(sum(test2(:,:),2));

for k=1:K
    term3(k)=sum(test3(k,(pilotIndex(k)==pilotIndex)'),2);
end

for k=1:K
    for s=1:tau_d
        SINR(k,s)=rho(s)^2*term1(k)/(term2(k)+rho(s)^2*term3(k)+term4(k));
    end
end

SE=squeeze(sum(real(log2(1+SINR(:,:))),2))/tau_c;
 
