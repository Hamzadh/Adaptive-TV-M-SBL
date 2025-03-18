function [X_mmse,D_mmse,theo,C] = MMSE_estimator(R,B,Phi,X,sig2e)
M=size(B,2);X_mmse=zeros(size(X));

for i=1:size(Phi,2)
    Y=B.'*conj(Phi(:,i));
    for l=1:size(Phi,2)
    Q_i(l,:,:)=abs(Phi(:,l).'*conj(Phi(:,i)))^2*squeeze(R(l,:,:));
    end
    Q(:,:)=inv(squeeze(sum(Q_i,1))+sig2e*eye(M));
    C(i,:,:)=squeeze(R(i,:,:))-squeeze(R(i,:,:))*Q*squeeze(R(i,:,:));
    X_mmse(i,:)=squeeze(R(i,:,:))*Q*Y;
    X_ls(i,:)=Y;
    theo_mmse(i)=1-trace(squeeze(R(i,:,:))*Q*squeeze(R(i,:,:)))/trace(squeeze(R(i,:,:)));
    %D_mmse_ue(i)= norm(X(i,:) - (X_mmse(i,:)))^2/norm(X(i,:))^2;
end
theo=mean(real(theo_mmse));

%Ru(:,:)=sum(R,1);P=pinv(Phi.');
%YY=B.'*P;
%Xmm=Ru*inv(Ru+sig2e*norm(P*P','fro')^2*eye(M))*YY;
% for i=1:size(Phi,2)
%     
%     Y=Phi(:,i)'*B;
%     for l=1:size(Phi,2)
%         
%         C_i(l,:,:)=abs(Phi(:,i)'*(Phi(:,l)))^2*(R(l,:,:));
%         
%     end
%     C(:,:)=sum(C_i,1);
%    % X_mmse(i,:)=squeeze(R(i,:,:))*inv(C+sig2e*eye(M))*Y.';
%      X_mmse(i,:)=squeeze(R(i,:,:))*(C+sig2e*eye(M))\Y.';
% 
%     %D_mmse= norm(X - (X_mmse),'fro')^2/norm(X,'fro')^2;
% 
% end

D_mmse= norm(X -X_mmse,'fro')^2/norm(X,'fro')^2;
%D_ls= norm(X -Xmm.','fro')^2/norm(X,'fro')^2

%DDD=mean(D_mmse_ue);

end

