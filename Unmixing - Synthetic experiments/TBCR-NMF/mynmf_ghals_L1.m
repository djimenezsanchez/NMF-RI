function [C,S,time,cost_record]= mynmf_ghals_L1(Y,C0,S0,phi)
%% HALSps+L1范数
%% =======================================================================
t0=clock;
Y(Y< 0) = eps;C=C0;S=S0;
J=size(C0,2);
tol=1e-4;
No_iter=1000;
%%
%lambda_sq=[];

cost_record=[];Flag_lambda=[];Sum_lambda=[];
inner=ones(1,J);%inner=[1 1 1];
Flag_sp=ones(1,J);Flag_sp(J)=0;%Flag_sp=[1 1 0];
lambda=0.01*ones(1,J);lambda(J)=0;%lambda=[0.01 0.01 0];
lambda_Fix=zeros(1,J);lambda_Fix(J)=1;%lambda_Fix=[0 0 1];
Count_small_incrase=0;Flag_Small_Increase=1;Small_increase=20;

cost = costfunction;
for k = 1: No_iter
    % 记录cost的变化
    cost_record=[cost_record,costfunction];
    % Update for A
    % 记录lambda的变化
    % lambda_sq=[lambda_sq;lambda]; figure(99);plot(lambda_sq);
    for j=J:-1:1%1:J
        C = nmfcore_L1(Y',S',C',j,Flag_sp(j),lambda(j))';
        if ~lambda_Fix(j)
            Flag_lambda=[Flag_lambda;zeros(1,J)];
            if Sparsity_Vec(C(:,j),0)<phi
                lambda(j)=lambda(j)*1.05;Flag_lambda(k,j)=-1;
            else
                if Flag_Small_Increase %% 令lambda导致的Sp略大于预期的phi，对于改善仿真数据中存在的局部性过强(AF存在空洞有帮助)
                    Count_small_incrase=Count_small_incrase+1;
                    if Count_small_incrase>Small_increase
                        Flag_Small_Increase=0;
                    end;
                    lambda(j)=lambda(j)*1.05;Flag_lambda(k,j)=-1;
                else
                    lambda(j)=lambda(j)*0.95;Flag_lambda(k,j)=1;
                end;
            end;
        end;
        S = nmfcore_L1(Y,C,S,j,0,0);
%         S=S./repmat(sqrt(sum(S.^2,2)),[1 size(S,2)]);
%         C=C.*repmat(sqrt(sum(S.^2,2))',[size(C,1) 1]);
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
        break;%% 当lambda全部陷于不稳定状态时，迭代终止
    end;
    
    checkstoppingcondition
    if stop ,
        break;
    end
    
end % k
time=etime(clock,t0);

    function X = nmfcore_L1(Y,A,X,j,Flag_sp,lambda)
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
        cost=sum(sum(Y-C*S))+lambda*sum(C)';
    end

%%   没有利用的部分
%     function X = nmfcore_phi(Y,A,X,j,Flag_sp,phi) % ,lcorr ,lsmth)
%         %Xsm = 0.001 * [X(:,2) (X(:,1:end-2)+X(:,3:end ))/2 X(:,end-1)];
%         AA = A'*A;
%         normA = diag (1./diag(AA));
%         AA = normA*AA;
%         AAX = normA*(A'*Y) - AA * X;%+ Xsm;
%         Xn= max(eps, X(j,:) + AAX(j,:));
%         if Flag_sp
%             Xn_max=max(Xn);
%             mybeta=fminsearch(@(mybeta) abs(Sparsity_Vec(Xn,mybeta)-phi),10);
%             Xn(Xn<(Xn_max/mybeta))=0;
%         end;
%         AAX = AAX - AA(:,j) * (Xn - X(j,:)) ;
%         Xn=max(Xn,eps);
%         X(j,:) = Xn;
%     end

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
