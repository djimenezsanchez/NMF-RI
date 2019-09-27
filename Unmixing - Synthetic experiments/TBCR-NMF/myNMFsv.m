function [C,S,t,objhistory]=myNMFsv(D,phi,C0,S0)%,tol)
t0=clock;
tol=1e-4;
maxiter=2000;   %%前两次数据  100
n=size(D,2);
% Data dimensions
vdim = size(D,1);
samples = size(D,2);
Comp=2;

C=C0;S=S0;

% Calculate initial objective
obj=sum(sum(D.^2));
objhistory = sum(sum((D-C*S).^2));

Flag_tol=1;
iter=0;
while iter<maxiter
    %figure(1);plot(S');figure(2);plot(C(:));
    % Update using standard NMF multiplicative update rule
    S = S.*(C'*D)./(C'*C*S + 1e-9);
    
    % Renormalize so rows of H have constant energy
    norms = sqrt(sum(S'.^2));
    S = S./(norms'*ones(1,samples));
    C = C.*(ones(vdim,1)*norms);
    
    C = (C).*(D*S')./(C*(S*S' )+ 1e-9);
    C_fluo=C(:,1:Comp);
    %% 计算Sparsity
                
    for inner=1:Comp
        C_fluo_Vec=C_fluo(:,inner);
        Sparsity_Vec(C_fluo_Vec,0);
        if Sparsity_Vec(C_fluo_Vec,0)<phi
            C_fluo_max=max(C_fluo_Vec);
            beta=fminsearch(@(beta) abs(Sparsity_Vec(C_fluo_Vec,beta)-phi),10);
            [Sparsity_Vec(C_fluo_Vec,0) C_fluo_max/beta];
            C_fluo_Vec(C_fluo_Vec<(C_fluo_max/beta))=0;
            C_fluo(:,inner)=C_fluo_Vec;
        end
    end;
    C(:,1:Comp)=C_fluo;
    
    % Update iteration count
    iter = iter+1;
    
    % Calculate objective
    newobj = sum(sum((D-C*S).^2))';
    objhistory = [objhistory newobj];
    
    if abs(objhistory(iter)-newobj)<tol
        %% 最后的处理
        S = S.*(C'*D)./(C'*C*S + 1e-9);
        % Renormalize so rows of H have constant energy
        norms = sqrt(sum(S'.^2));
        S = S./(norms'*ones(1,samples));
        C = C.*(ones(vdim,1)*norms);
        C = (C).*(D*S')./(C*(S*S' )+ 1e-9);
        %break;
    end;
    
end
t=etime(clock,t0);

    function Sparsity_Ap=Sparsity_Vec(Ap_Vec,beta)
        Ap_Dim=max(size(Ap_Vec,1),size(Ap_Vec,2));
        Ap_max=max(Ap_Vec);
        if (beta==0)
            threshold=0;
        else
            threshold=Ap_max/beta;
        end;
%         Ap_Vec(Ap_Vec<threshold)=1e-9;
        Ap_Vec(Ap_Vec<threshold)=0;
        %Ap_Vec(Ap_Vec>threshold)=Ap_Vec(Ap_Vec>threshold)-threshold;
        %Ap_Vec(Ap_Vec<0)=0;
        Ap_norm1=sum(abs(Ap_Vec));
        Ap_norm2=sqrt(sum(Ap_Vec.^2));
        Sparsity_Ap=(sqrt(Ap_Dim)-Ap_norm1/Ap_norm2)/(sqrt(Ap_Dim)-1);
    end
end