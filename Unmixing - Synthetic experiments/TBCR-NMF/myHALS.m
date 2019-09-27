function [C,S]= myHALS12(D,C0,S0)
maxIter = 1000;
eps = 1e-9;
[M,N] = size(C0);

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

        C(:,k) = max(eps, Rk*SS(k,:)'/(SS(k,:)*SS(k,:)'));
        S(k,:) = max(eps, C(:,k)'*Rk/(C(:,k)'*C(:,k)));

        S=S./repmat(sqrt(sum(S.^2,2)),[1 size(S,2)]);
        C=C.*repmat(sqrt(sum(S.^2,2))',[size(C,1) 1]);
        
    end
    CC = C;SS=S;
%     end
end


