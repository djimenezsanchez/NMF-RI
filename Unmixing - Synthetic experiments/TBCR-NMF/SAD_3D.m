%% SAD_3D
function [SAD_all]=SAD_3D(Pn,S_all)
%% 计算SAD，Pn为标准光谱，维度：组分数*Bands，S_all求得的光谱，第三维度表示次数或者不同参数下的值
[VD,~,Num]=size(S_all);
SAD_all=zeros(Num,VD);
for i=1:Num
    S=S_all(:,:,i);
    SAD_all(i,:)=SAD(Pn,S);
end

    function [ SAD_Value ] = SAD( S1,S2 )
        %UNTITLED2 Summary of this function goes here
        %   Detailed explanation goes here
        %  spectral angle measure/Distance   SAD/SAD
        %  S1 S2: Component Number*Spectral band
        %pixels=size(S1,1);
        S1_norm2=sqrt(sum(S1.^2,2));
        S2_norm2=sqrt(sum(S2.^2,2));
        SAD_Value=real(acos((sum(S1.*S2,2))./(S1_norm2.*S2_norm2)));
        SAD_Value=SAD_Value';
    end
end