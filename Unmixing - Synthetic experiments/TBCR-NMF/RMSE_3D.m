function [RMSE_all]=RMSE_3D(An,C_all)
[~,VD,Num]=size(C_all);
RMSE_all=zeros(Num,VD);
for i=1:Num
    C=C_all(:,:,i);
    RMSE_all(i,:)=RMSE(An,C);
end

    function [RMSE_Value] = RMSE(C1,C2)
        %UNTITLED5 Summary of this function goes here
        %   Detailed explanation goes here
        %  C1 C2: Pixels * Component Number
        component=3;
        for j=1:component
            RMSE_Value(1,j)=sqrt(mean((C1(:,j)-C2(:,j)).^2));
        end;
        %RMSE_Value(1,4)=mean(RMSE_Value(1,[1 2 3]));
    end
end