function [eta_opt] = maxmin_power_control(SE_MF,R,gamma_kl,K,L,tau_c,tau_d,Pmax,pilotIndex,rho)

tmin=2^(min(SE_MF(:)))-1;
tmax=2^(2*min(SE_MF)+1.2)-1;
epsi=max(tmin/5,0.01);




sigma2=1;
p=Pmax;

a_kl=ones(K,L);

b_k=zeros(K,L);
Y=zeros(L,L,K,K);
c_ki=zeros(K,K,L);
A=zeros(L,L,K);

test1=zeros(1,K);
test2=zeros(K,K);
test3=zeros(K,K);
test4=zeros(1,K);


a_k=a_kl(:,:);

for l=1:L
    for k=1:K
        b_k(k,l)=real(trace(gamma_kl(:,:,k,l)));
        A(l,l,k)=real(trace(gamma_kl(:,:,k,l)));
        for i=1:K
            Y(l,l,k,i)=real(trace(gamma_kl(:,:,k,l)*R(:,:,i,l)));
            c_ki(k,i,l)=sqrt(real(trace(gamma_kl(:,:,k,l)*gamma_kl(:,:,i,l))));
        end
    end
end

for k=1:K
    test1(k)=abs(a_k(k,:)*b_k(k,:)')^2;
    test4(k)=sigma2*a_k(k,:)*A(:,:,k)*a_k(k,:)';
    for i=1:K
        test2(k,i)=a_k(k,:)*Y(:,:,k,i)*a_k(k,:)';
        if i==k
           test3(k,i)=0;
        else
           test3(k,i)=abs(a_k(k,:)*squeeze(c_ki(k,i,:)))^2; 
        end
    end
end


%cvx_solver sedumi
 cvx_quiet true
            while( tmax - tmin > epsi)

            tnext = (tmax+tmin)/2; 
           cvx_begin %sdp
              variables x(K,1) 
              minimize(0)
              subject to
                for k=1:K
                    p*test2(k,:)*x(:)+p*(test3(k,:).*(pilotIndex(k)==pilotIndex))*x(:)+test4(k) <= (1/tnext)*p*x(k)*test1(k);                                 
                end              
                for k=1:K
                    x(k)>=0;
                    x(k)<=1;
                end
                            
            cvx_end

            % bisection
            if contains(cvx_status,'Solved') % feasible
            disp('Problem is feasible'); 
            %fprintf(1,'Problem is feasible ',tnext);
            tmin = tnext;
            eta_opt=x;
            else % not feasible
            disp('Problem is not feasible'); 
            %fprintf(1,'Problem not feasible ',tnext);
            tmax = tnext;   
            end

            end




