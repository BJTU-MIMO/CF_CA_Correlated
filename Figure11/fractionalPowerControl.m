function [p_CF_control,p_SC_control] = fractionalPowerControl(beta,K,L,P,theta_CF,theta_SC)

p_CF_control=zeros(K,L);
p_SC_control=zeros(K,L);
index=zeros(1,K);
p_SC_k=zeros(1,K);

p_CF_k=squeeze(1./(sum(beta(:,:),2).^theta_CF));
p_CF_1=p_CF_k./(sum(p_CF_k(:)))*(K*P);
for k=1:K
    for l=1:L
        p_CF_control(k,l)=p_CF_1(k);
    end
end


for k=1:K
[~,index(k)]= max(beta(k,:));   
end

for k=1:K
    p_SC_k(k)=1/(beta(k,index(k))^theta_SC);
end
p_SC_1=p_SC_k(:)./(sum(p_SC_k(:)))*(K*P);
for k=1:K
    for l=1:L
        p_SC_control(k,l)=p_SC_1(k);
    end
end