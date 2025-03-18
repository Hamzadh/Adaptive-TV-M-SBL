function [SBL_xhat,t,mse_SBL,SBL_SRR,warm_init] = sbl_mmv(A,y,N,true_x,true_support,sig2e,Tau_p)
max_iter_sbl=351;M=size(y,2);

stopping_creterion=1e-2;
%Sigma_y=A*diag(gamma*A';
    Gamma(:,:,1)=.1*eye(N);
SBL_SRR=zeros(1);
for t=1:max_iter_sbl
    Sigma_y=(A*Gamma(:,:,t)*A'+sig2e*eye(Tau_p));
    mu_x=Gamma(:,:,t)*A'*inv(Sigma_y)*y;
   Sigma_x=Gamma(:,:,t)-Gamma(:,:,t)*A'*inv(Sigma_y)*A*Gamma(:,:,t);
    for l=1:N
        gamma(l,t)=(1/M)*(norm(mu_x(l,:))^2)+Sigma_x(l,l);
    end
    Gamma(:,:,t+1)=diag(gamma(:,t));
    xhat_sbl(:,:,t)=mu_x;
    if t>150
    if norm(xhat_sbl(:,t)-xhat_sbl(:,:,t-1))<stopping_creterion
       break; 
    end
    end
end
SBL_xhat=xhat_sbl(:,:,t);
warm_init=gamma(:,min(t,200));%
mse_SBL=norm(SBL_xhat(true_support,:)-true_x(true_support,:),'fro')^2/norm(true_x(true_support,:),'fro')^2;
%mse_SBL=norm(SBL_xhat-true_x,'fro')^2/norm(true_x,'fro')^2;

support_sbl=find(abs(vecnorm(SBL_xhat,2,2))>=0.2)';

SBL_SRR=F1_score(SBL_xhat,true_support,.1);

% if isempty(setdiff(support_sbl,true_support))
%     SBL_SRR=1;
% else
%     SBL_SRR=0;
% end
%stem(real(gamma(:,t)))
end

