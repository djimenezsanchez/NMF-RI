function [A,S,conv]=hc_sparseNMF(V,phi,A0,S0,targetNum)%,tol)

tol=1e-4;
N_iter=1000; 
A = A0; S = S0;
[m, n] = size(V);


iter = 1;
conv = ones(1,N_iter);
h = waitbar(0,'开始迭代');
while iter < N_iter
    oldconv = sum(sum((V - A*S).^2));
    
    S = S.*(A'*V)./((A'*A)*S + 1e-9); 
    % Renormalize so rows of S have constant energy
    norms = sqrt(sum(S.^2,2));
    S = S./repmat(norms, 1, n);
    A = A.*(ones(m,1)*norms');
    A = A.*(V*S')./(A*(S*S') + 1e-9);
    
    
    for i = 1:targetNum
        updateA = A(:,i);
        if Sparsity_Vec(updateA,0) < phi
%         beta=fminsearch(@(beta) abs(Sparsity_Vec(updateA,beta)-phi),10,optimset('MaxIter',100,'TolX',0.01));
            beta=fminsearch(@(beta) abs(Sparsity_Vec(updateA,beta)-phi),100);
            beta
            updateA(updateA < max(updateA)/beta) = 1e-9;
        end
        A(:,i) = updateA;
    end
    
    newconv = sum(sum((V - A*S).^2));
    conv(1,iter) = newconv;
    if abs(oldconv - newconv) < tol
        break;
    end
    
    Row = 256;Col = 256;
%     figure(100),imagesc(reshape(A(:,1),Row,Col));
%     figure(101),imagesc(reshape(A(:,2),Row,Col));
%     figure(102),imagesc(reshape(A(:,3),Row,Col));
%     figure(103),plot(conv);
    figure(100),imagesc([reshape(A(:,1),Row,Col) reshape(A(:,2),Row,Col) reshape(A(:,3),Row,Col)]);
    figure(101),plot(S');
    iter = iter + 1;
    
    waitbar(iter/1000,h,'正在迭代...');
end

close(h);

S = S.*(A'*V)./((A'*A)*S + 1e-9); 
% Renormalize so rows of S have constant energy
norms = sqrt(sum(S.^2,2));
S = S./repmat(norms, 1, n);
A = A.*(ones(m,1)*norms');

A = A.*(V*S')./(A*(S*S') + 1e-9);


end

