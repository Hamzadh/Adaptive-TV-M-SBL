function [SBL_xhat,t,mse_SBL,SBL_SRR,warm_init] = pattern_coupled_MSBL(A,y,N,true_x,true_support,sig2e,Tau_p,init)
max_iter_sbl=201;M=size(y,2);M=1;
stopping_creterion=1e-4;
%Sigma_y=A*diag(gamma*A';
Gamma(:,:,1)=(diag(init)*eye(N));gamma_inv(:,1)=diag(inv(Gamma(:,:,1)));SBL_SRR=zeros(1);
a=1.5;b=0.00001;beta=1;
for t=1:max_iter_sbl
    Sigma_y=(A*Gamma(:,:,t)*A'+sig2e*eye(Tau_p));
    mu_x=Gamma(:,:,t)*A'/(Sigma_y)*y;
    Sigma_x=Gamma(:,:,t)-Gamma(:,:,t)*A'/(Sigma_y)*A*Gamma(:,:,t);
   % Sigma_x=zeros(N);
    for l=1:N
        if l==1
            w(l)= real(norm(mu_x(l,:))^2/M+Sigma_x(l,l))+beta*real((norm(mu_x(l+1,:))^2/M)+Sigma_x(l+1,l+1));
            gamma_inv(l,t)=real((a)/(.5*w(l)+b));
        elseif l==N
            w(l)=real( ((norm(mu_x(l,:))^2/M)+Sigma_x(l,l))+beta*real(((norm(mu_x(l-1,:))^2/M)+Sigma_x(l-1,l-1))));
            gamma_inv(l,t)=real((a)/(.5*w(l)+b));
        else
            w(l)= real(norm(mu_x(l,:))^2/M+Sigma_x(l,l))+beta*real(norm(mu_x(l+1,:))^2/M+Sigma_x(l+1,l+1))+...
                +beta*real((norm(mu_x(l-1,:))^2/M)+Sigma_x(l-1,l-1)) ;
            gamma_inv(l,t)=real((a)/(.5*w(l)+b));
        end
        if gamma_inv(l,t)>10^7
            gamma_inv(l,t)=10^7;
        end        
    end
    for i=1:N
        if i==1
               Gamma(i,i,t+1)=   gamma_inv(i,t)+beta*gamma_inv(i+1,t);

        elseif i==N
         Gamma(i,i,t+1)=   gamma_inv(i,t)+beta*gamma_inv(i-1,t);
        else
       Gamma(i,i,t+1)=   gamma_inv(i,t)+beta*gamma_inv(i-1,t)+beta*gamma_inv(i+1,t);
        end
    end
    Gamma(:,:,t+1)=inv(Gamma(:,:,t+1));
    xhat_sbl(:,:,t)=mu_x;
    if t>20
        if norm(xhat_sbl(:,:,t)-xhat_sbl(:,:,t-1))<stopping_creterion
            break;
        end
    end
    mse_SdBL(t)=norm(xhat_sbl(:,:,t)-true_x,'fro')^2/norm(true_x,'fro')^2;
end
SBL_xhat=xhat_sbl(:,:,t);
warm_init=xhat_sbl(:,:,10);
mse_SBL=norm(SBL_xhat(true_support,:)-true_x(true_support,:),'fro')^2/norm(true_x(true_support,:),'fro')^2;
%
support_sbl=find(abs(vecnorm(SBL_xhat,2,2))>=0.2)';
SBL_SRR=F1_score(SBL_xhat,true_support,.2);

% if isempty(setdiff(support_sbl,true_support))
%     SBL_SRR=1;
% else
%     SBL_SRR=0;
% end
%stem(real(gamma(:,t)))
end



