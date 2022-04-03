function [beta,h,hhat,hc,gn,gamma_kl,Np] = functionChannelEstimates(R,gain,K,L,N,nbrOfRealizations,p,pilotIndex,rho)

beta=gain(:,:);

%% Generate the channel realizations
h=zeros(N,nbrOfRealizations,K,L);
gn=zeros(N,nbrOfRealizations,K,L);

W = sqrt(0.5)*(randn(N,nbrOfRealizations,K,L)+1i*randn(N,nbrOfRealizations,K,L));
Wn= sqrt(0.5)*(randn(N,nbrOfRealizations,K,L)+1i*randn(N,nbrOfRealizations,K,L));

for k=1:K
    for l=1:L
        Rsqrt = sqrtm(R(:,:,k,l));
        h(:,:,k,l)=Rsqrt*W(:,:,k,l);
        gn(:,:,k,l)=Rsqrt*Wn(:,:,k,l);
    end
end

%% Perform channel estimation
sigma2=1;
eyeN = eye(N);
Psi=zeros(N,N,K,L);
gamma_kl=zeros(N,N,K,L);
%Generate realizations of normalized noise
Np = sqrt(0.5*sigma2)*(randn(N,nbrOfRealizations,K,L) + 1i*randn(N,nbrOfRealizations,K,L));

hhat=zeros(N,nbrOfRealizations,K,L);
z=zeros(N,nbrOfRealizations,K,L);


for k=1:K
    for l=1:L
        %z(:,:,k,l)=sqrt(p)*sum(h(:,:,find(pilotIndex(k)==pilotIndex)',l),3)+Np(:,:,k,l);
        z(:,:,k,l)=sqrt(p)*rho(pilotIndex(k))*sum(h(:,:,find(pilotIndex(k)==pilotIndex)',l),3)+sqrt(p)*sqrt(1-rho(pilotIndex(k))^2)*sum(gn(:,:,find(pilotIndex(k)==pilotIndex)',l),3)+Np(:,:,k,l);
        Psi(:,:,k,l)=p*sum(R(:,:,find(pilotIndex(k)==pilotIndex)',l),3)+sigma2*eyeN;
        hhat(:,:,k,l)=rho(pilotIndex(k))*sqrt(p)*(R(:,:,k,l)/Psi(:,:,k,l))*z(:,:,k,l);
        gamma_kl(:,:,k,l)=rho(pilotIndex(k))^2*p*(R(:,:,k,l)/Psi(:,:,k,l))*R(:,:,k,l);
    end
end

hc=h(:,:,:,:)-hhat(:,:,:,:);




