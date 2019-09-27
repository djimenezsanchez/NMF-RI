function [C,S,time,cost_record,lambda_sq,k]= mynmf_ghals_L12_Ultimate(Y,C0,S0,phi)
%% HALSps+L0.5 norm
%% =======================================================================
t0=clock;
Ratio=0.8;
C0=C0+0.1*repmat(max(C0),[size(C0,1) 1]);
 %C0(:,1:2)=Ratio*C0(:,1:2);
Y(Y< 0) = eps;C=C0;S=S0;
J=size(C0,2);
tol=1e-3;
No_iter=500;
%%
lambda_sq=[];

cost_record=[];Flag_lambda=[];Sum_lambda=[];
inner=ones(1,J);%inner=[1 1 1];
Flag_sp=ones(1,J);Flag_sp(J)=0;%Flag_sp=[1 1 0];
lambda=0.01*ones(1,J);lambda(J)=0;%lambda=[0.01 0.01 0];
lambda_Fix=zeros(1,J);lambda_Fix(J)=1;%lambda_Fix=[0 0 1];
Count_small_incrase=0;Flag_Small_Increase=1;Small_increase=0;

%phi_true=phi;
%phi=0.85;


cost = costfunction;Once=1;Second=1;
for k = 1: No_iter
    % record the change of cost
    cost_record=[cost_record,costfunction];
    % Update for A
    % record the change of lambda
      lambda_sq=[lambda_sq;lambda]; %figure(99);plot(lambda_sq);
     % figure(1);plot(S');figure(2);plot(C(:))
    


        for j=J:-1:1%% previous version 1:J 4_24
            C = nmfcore_L12(Y',S',C',j,Flag_sp(j),lambda(j))';
            if ~lambda_Fix(j)
                Flag_lambda=[Flag_lambda;zeros(1,J)];
                if Sparsity_Vec(C(:,j),0)<phi
                    lambda(j)=lambda(j)*1.05;Flag_lambda(k,j)=-1;
                else
                    if Flag_Small_Increase %% make the lambda-resulted sparsity be greater than the expectant phi£¬it is helpful for the hole of background AF
                        Count_small_incrase=Count_small_incrase+1;
                        if Count_small_incrase>Small_increase
                            Flag_Small_Increase=0;
                        end;
                        lambda(j)=lambda(j)*1.1;Flag_lambda(k,j)=-1;
                    else
                        %phi=phi_true;
                        lambda(j)=lambda(j)*0.95;Flag_lambda(k,j)=1;
%                         if lambda(j)<1e-12
%                             lambda_Fix(j)=1
%                         end;
                    end;
                end;
            end;
        end;
       
        Size=128;
        %figure(22);imagesc(reshape(C(:,3),[Size Size]));axis image
       % figure(33);imagesc(reshape(C(:,2),[Size Size]));axis image
        %figure(11);imagesc(reshape(C(:,1),[Size Size]));axis image

        if sum(Flag_lambda(k,1:2))==2 && Once
            C(:,1:2)=Ratio*C(:,1:2);
            Once=0
        end;
        
        for j=J:-1:1
            S = nmfcore_L1(Y,C,S,j,0,0);
        end;
    S=S./repmat(sqrt(sum(S.^2,2)),[1 size(S,2)]);
    C=C.*repmat(sqrt(sum(S.^2,2))',[size(C,1) 1]);
  

    if max(inner)+1>size(Sum_lambda,1)
        Sum_lambda=[Sum_lambda;zeros(1,J)];
    end;
    for j=1:J
        if ~lambda_Fix(j) && k>0 && mod(k,20)==0
            Sum_lambda(inner(j),j)=sum(Flag_lambda(20*(inner(j)-1)+1:20*inner(j),j),1);
            if inner(j)>3
                temp=Sum_lambda(inner(j)-2:inner(j),j);
                if abs(temp(1))~=20 && abs(temp(2))~=20
                    if abs(temp(1)-temp(2))<5 && abs(temp(2)-temp(3))<5
                        lambda_Fix(j)=1;
                    end;
                end;
            end;
            inner(j)=inner(j)+1;
        end;
    end;
    if sum(lambda_Fix(1:J-1))~=J-1
        lambda_Fix(1:J-1)=0;
    else
        break;%% when lambda is in the unstable state£¬the iteration stops
    end;
    
    checkstoppingcondition
    if stop ,
        break;
    end
    
end % k

% iter=1;maxiter=500;Comp=2;D=Y;objhistory=sum(sum((D-C*S).^2))';S(3,:)=S0(3,:);
% while iter<maxiter
% %     for j=J:-1:1%% previous version 1:J 4_24
% %         C = nmfcore_phi(Y',S',C',j,Flag_sp(j),phi)';
% %         C=max(C,0);
% %         S = nmfcore_phi(Y,C,S,j,0,0);
% %         %         S=S./repmat(sqrt(sum(S.^2,2)),[1 size(S,2)]);
% %         %         C=C.*repmat(sqrt(sum(S.^2,2))',[size(C,1) 1]);
% %     end;
% %     S=S./repmat(sqrt(sum(S.^2,2)),[1 size(S,2)]);
% %     C=C.*repmat(sqrt(sum(S.^2,2))',[size(C,1) 1]);
%         figure(1);plot(S')
%         % Update using standard NMF multiplicative update rule
%
%         S = S.*(C'*D)./(C'*C*S + 1e-9);
%         %if iter==1
%             S(3,:)=S0(3,:);
%         %end
%         %S=S0;
%
%         % Renormalize so rows of H have constant energy
%         S=S./repmat(sqrt(sum(S.^2,2)),[1 size(S,2)]);
%         C=C.*repmat(sqrt(sum(S.^2,2))',[size(C,1) 1]);
%
%         C = (C).*(D*S')./(C*(S*S' )+ 1e-9);
%         C_fluo=C(:,1:Comp);
%         %% compute the sparsity
%         for inner=1:Comp
%             C_fluo_Vec=C_fluo(:,inner);
%             Sparsity_Vec(C_fluo_Vec,0)
%             %if Sparsity_Vec(C_fluo_Vec,0)<phi
%                 C_fluo_max=max(C_fluo_Vec);
%                 beta=fminsearch(@(beta) abs(Sparsity_Vec(C_fluo_Vec,beta)-phi),10);
%                 C_fluo_Vec(C_fluo_Vec<(C_fluo_max/beta))=0;
%                 C_fluo(:,inner)=C_fluo_Vec;
%             %end
%         end;
%         C(:,1:Comp)=C_fluo;
%
%     % Update iteration count
%     iter = iter+1;
%
%     % Calculate objective
%     newobj = sum(sum((D-C*S).^2))';
%     objhistory = [objhistory newobj];
%
% %     if abs(objhistory(iter)-newobj)<tol
% %         %% the last processing
% %         S = S.*(C'*D)./(C'*C*S + 1e-9);
% %         % Renormalize so rows of H have constant energy
% %         S=S./repmat(sqrt(sum(S.^2,2)),[1 size(S,2)]);
% %         C=C.*repmat(sqrt(sum(S.^2,2))',[size(C,1) 1]);
% %         C = (C).*(D*S')./(C*(S*S' )+ 1e-9);
% %         break;
% %     end;
%
% end


time=etime(clock,t0);

    function X = nmfcore_L1(Y,A,X,j,Flag_sp,lambda)
        %% L1 norm
        AA = A'*A;
        normA = diag (1./diag(AA));
        AA = normA*AA;
        AAX = normA*(A'*Y) - AA * X;
        if Flag_sp
            
            Xn= max(eps, X(j,:) + AAX(j,:)-lambda);
        else
            Xn= max(eps, X(j,:) + AAX(j,:));
        end;
        AAX = AAX - AA(:,j) * (Xn - X(j,:)) ;
        X(j,:) = Xn;
    end

    function X = nmfcore_L12(Y,A,X,j,Flag_sp,lambda)
        %% L0.5 norm
        AA = A'*A;
        normA = diag (1./diag(AA));
        AA = normA*AA;
        AAX = normA*(A'*Y) - AA * X;
        if Flag_sp
            Xn= max(eps, X(j,:) + AAX(j,:)-lambda.*power(X(j,:)+eps,-0.5));
        else
            Xn= max(eps, X(j,:) + AAX(j,:));
        end;
        AAX = AAX - AA(:,j) * (Xn - X(j,:)) ;
        X(j,:) = Xn;
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

    function checkstoppingcondition
        cost_old = cost; cost = costfunction;
        %stop = abs(cost_old-cost) <= tol*cost;
        stop = abs(cost_old-cost) <= tol;
    end

    function cost = costfunction
        %Yhat = A*X+eps;
        %cost = sum(sum(Y.*log(Y./Yhat + eps) - Y + Yhat));
        cost=sum(sum(Y-C*S))+lambda*sum(power(C+eps,0.5))';
    end

%%   the deprecated part of algorithm
    function X = nmfcore_phi(Y,A,X,j,Flag_sp,phi) % ,lcorr ,lsmth)
        %Xsm = 0.001 * [X(:,2) (X(:,1:end-2)+X(:,3:end ))/2 X(:,end-1)];
        AA = A'*A;
        normA = diag (1./diag(AA));
        AA = normA*AA;
        AAX = normA*(A'*Y) - AA * X;%+ Xsm;
        Xn= max(eps, X(j,:) + AAX(j,:));
        if Flag_sp
            Xn_max=max(Xn);
            mybeta=fminsearch(@(mybeta) abs(Sparsity_Vec(Xn,mybeta)-phi),10);
            Xn(Xn<(Xn_max/mybeta))=0;
        end;
        AAX = AAX - AA(:,j) * (Xn - X(j,:)) ;
        Xn=max(Xn,0);
        X(j,:) = Xn;
    end

%     function Sparsity_Ap=Sparsity_Vec(Ap_Vec,beta)
%         Ap_Dim=max(size(Ap_Vec,1),size(Ap_Vec,2));
%         Ap_max=max(Ap_Vec);
%         if (beta==0)
%             threshold=0;
%         else
%             threshold=Ap_max/beta;
%         end;
%         Ap_Vec(Ap_Vec<threshold)=0;
%         %Ap_Vec(Ap_Vec>threshold)=Ap_Vec(Ap_Vec>threshold)-threshold;
%         %Ap_Vec(Ap_Vec<0)=0;
%         Ap_norm1=sum(abs(Ap_Vec));
%         Ap_norm2=sqrt(sum(Ap_Vec.^2));
%         Sparsity_Ap=(sqrt(Ap_Dim)-Ap_norm1/Ap_norm2)/(sqrt(Ap_Dim)-1);
%     end
end
%end
