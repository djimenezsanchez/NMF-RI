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