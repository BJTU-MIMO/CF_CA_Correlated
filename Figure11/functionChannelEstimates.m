function [beta,h,hhat,hc,gn,gamma_kl,Np] = functionChannelEstimates(gain,K,L,nbrOfRealizations,p,pilotIndex,rho,rhohat)

beta=gain(:,:);

%% Generate the channel realizations
h=zeros(nbrOfRealizations,K,L);
gn=zeros(nbrOfRealizations,K,L);


W = (randn(nbrOfRealizations,K,L)+1i*randn(nbrOfRealizations,K,L));
Wn= (randn(nbrOfRealizations,K,L)+1i*randn(nbrOfRealizations,K,L));

for k=1:K
    for l=1:L
        h(:,k,l)=sqrt(0.5*beta(k,l))*W(:,k,l);
        gn(:,k,l)=sqrt(0.5*beta(k,l))*Wn(:,k,l);
    end
end

%% Perform channel estimation
N0=1;
ss=zeros(K,L);
gamma_kl=zeros(K,L);
%Generate realizations of normalized noise
Np = sqrt(0.5*N0)*(randn(nbrOfRealizations,K,L) + 1i*randn(nbrOfRealizations,K,L));

hhat=zeros(nbrOfRealizations,K,L);
z=zeros(nbrOfRealizations,K,L);


for k=1:K
    for l=1:L
        z(:,k,l)=sqrt(p)*sum(rho(pilotIndex(k))*h(:,find(pilotIndex(k)==pilotIndex)',l),2)+sqrt(p)*sum(rhohat(pilotIndex(k))*gn(:,find(pilotIndex(k)==pilotIndex)',l),2)+Np(:,k,l);
        ss(k,l)=rho(pilotIndex(k))*sqrt(p)*beta(k,l)/(p*sum(beta(find(pilotIndex(k)==pilotIndex)',l),1)+N0);
        hhat(:,k,l)=ss(k,l)*z(:,k,l); 
        gamma_kl(k,l)=rho(pilotIndex(k))^2*p*beta(k,l)^2/(p*sum(beta(find(pilotIndex(k)==pilotIndex)',l),1)+N0);
    end
end

hc=h(:,:,:)-hhat(:,:,:);




