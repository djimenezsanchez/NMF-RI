function [C,S]= myHALS12(D,C0,S0,phi)
maxIter = 1000;
eps = 1e-9;
[M,N] = size(C0);
lambda = 0.001*ones(1,N-1); 
cost = [];


CC = C0;SS=S0;
C = C0;S=S0;

for iter = 1:maxIter
    for k = 1:N
        % Rk
        Rk = D;
        for j = 1:N
            if j ~= k
                Rk = Rk - CC(:,j)*SS(j,:);
            end
        end
        % target
        if k ~= N
            C(:,k) = max(eps, (Rk*SS(k,:)' - lambda(1,k)*(1./sqrt(CC(:,k))))/(SS(k,:)*SS(k,:)'));
            if Sparsity_Vec(C(:,k),0)> phi + 0.001
                lambda(1,k) = 0.95*lambda(1,k);
            elseif Sparsity_Vec(C(:,k),0)< phi - 0.001
                lambda(1,k) = 1.05*lambda(1,k);
            else
            end
        else
            C(:,k) = max(eps, Rk*SS(k,:)'/(SS(k,:)*SS(k,:)'));
        end
        S(k,:) = max(eps, C(:,k)'*Rk/(C(:,k)'*C(:,k)));

        S=S./repmat(sqrt(sum(S.^2,2)),[1 size(S,2)]);
        C=C.*repmat(sqrt(sum(S.^2,2))',[size(C,1) 1]);
        
    end
%     newnorm = sum(sum((D - C*S).^2));
%     oldnorm = sum(sum((D - CC*SS).^2));cost = [cost  oldnorm];
%     if abs(newnorm - oldnorm) < 1e-4
%         break;
%     else
        CC = C;SS=S;
%     end
    
    figure(100),imagesc([reshape(C(:,1), sqrt(M), sqrt(M)),reshape(C(:,2), sqrt(M), sqrt(M)),reshape(C(:,3), sqrt(M), sqrt(M))] );
    figure(101),plot(cost);
end

     function Sparsity_Ap=Sparsity_Vec(Ap_Vec,beta)
        Ap_Dim=max(size(Ap_Vec,1),size(Ap_Vec,2));
        Ap_max=max(Ap_Vec);
        if (beta==0)
            threshold=0;
        else
            threshold=Ap_max/beta;
        end;
        Ap_Vec(Ap_Vec<threshold)=0;
        %Ap_Vec(Ap_Vec>threshold)=Ap_Vec(Ap_Vec>threshold)-threshold;
        %Ap_Vec(Ap_Vec<0)=0;
        Ap_norm1=sum(abs(Ap_Vec));
        Ap_norm2=sqrt(sum(Ap_Vec.^2));
        Sparsity_Ap=(sqrt(Ap_Dim)-Ap_norm1/Ap_norm2)/(sqrt(Ap_Dim)-1);
     end
end

