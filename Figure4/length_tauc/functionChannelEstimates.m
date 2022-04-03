function [gamma_kl,gamma_kil,tracQ] = functionChannelEstimates(R,K,L,N,p,pilotIndex,rho)

sigma2=1;
eyeN = eye(N);
Psi=zeros(N,N,K,L);
gamma_kl=zeros(N,N,K,L);
gamma_kil=zeros(N,N,K,K,L);
tracQ=zeros(K,L);

for k=1:K
    for l=1:L
        Psi(:,:,k,l)=p*sum(R(:,:,find(pilotIndex(k)==pilotIndex)',l),3)+sigma2*eyeN;
        gamma_kl(:,:,k,l)=rho(pilotIndex(k))^2*p*(R(:,:,k,l)/Psi(:,:,k,l))*R(:,:,k,l);
        tracQ(k,l)=abs(trace(gamma_kl(:,:,k,l)));
        for i=1:K
            gamma_kil(:,:,k,i,l)=rho(pilotIndex(k))*rho(pilotIndex(i))*p*(R(:,:,i,l)/Psi(:,:,k,l))*R(:,:,k,l);
        end
    end
end






