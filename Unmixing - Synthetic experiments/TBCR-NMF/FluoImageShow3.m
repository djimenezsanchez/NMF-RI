% Data为图像的原始数据
% BBmp为背景图片
% bgV为设置的背景阈值
% maxV为映射范围 映射范围为bgV~maxV
% minV < bV < maxV
function Imgshow=FluoImageShow3(Data,BBmp,bgV,maxV,a)
% Data=medfilt2(Data);
[m,n,chl]=size(BBmp);
Imgshow=zeros(m,n,chl);
minV=min(Data(:));

if bgV<minV
    error('The background value should be bigger than the minimum of the imported data ');
end


Data=(Data-bgV)/(maxV-bgV);% 归一化
Data(Data<0)=0;             % 小于背景值的部分为0
Data=log(Data*(a-1)+1)/log(a)*255;% 映射到0~255
% Data=medfilt2(Data);
r=1;b=2;g=3;c=255;
% figure(1),imshow(uint8(Data));

for i=1:m
    for j=1:n
        if Data(i,j)==0
            Imgshow(i,j,r)=BBmp(i,j,r);
            Imgshow(i,j,b)=BBmp(i,j,g);
            Imgshow(i,j,g)=BBmp(i,j,b);
        else
            if  Data(i,j)<=c/4
                Imgshow(i,j,r)=0;
                Imgshow(i,j,b)=4*Data(i,j);
                Imgshow(i,j,g)=c;
            elseif Data(i,j)<=c/2
                Imgshow(i,j,r)=0;
                Imgshow(i,j,b)=c;
                Imgshow(i,j,g)=-4*Data(i,j)+2*c;
            elseif Data(i,j)<=3*c/4
                Imgshow(i,j,r)=4*Data(i,j)-2*c;
                Imgshow(i,j,b)=c;
                Imgshow(i,j,g)=0;
            else
                Imgshow(i,j,r)=c;
                Imgshow(i,j,b)=-4*Data(i,j)+4*c;
                Imgshow(i,j,g)=0;
            end
        end
    end
end

                

% % 进行伪彩色变换
% P1=255/4;P2=255/2;P3=255*3/4;
% 
% DataCOLR=DataCOL;
% DataCOLR(DataCOLR<P2)=0;
% DataCOLR((DataCOLR>=P2)&(DataCOLR<=P3))=4*DataCOLR((DataCOLR>=P2)&(DataCOLR<=P3))-255*2;
% DataCOLR(DataCOLR>P3)=255;
% 
% DataCOLG=DataCOL;
% DataCOLG(DataCOLG<P1)=4*DataCOLG(DataCOLG<P1);
% DataCOLG((DataCOLG>=P1)&(DataCOLG<=P3))=255;
% DataCOLG(DataCOLG>P3)=-4*DataCOLG(DataCOLG>P3)+4*255;
% 
% DataCOLB=DataCOL;
% DataCOLB(DataCOLB<P1)=255;
% DataCOLB((DataCOLB>=P1)&(DataCOLB<=P2))=-4*DataCOLB((DataCOLB>=P1)&(DataCOLB<=P2))+2*255;
% DataCOLB(DataCOLB>P2)=0;
% 
% % 将变换后的3通道数据整合到一起
% I1(1:m,1:n)=0;
% I1(indexBG)=BackData(indexBG);
% I1(indexCOL)=DataCOLR;
% 
% I2(1:m,1:n)=0;
% I2(indexBG)=BackData(indexBG);
% I2(indexCOL)=DataCOLB;
% 
% 
% I3(1:m,1:n)=0;
% I3(indexBG)=BackData(indexBG);
% I3(indexCOL)=DataCOLG;
% 
% I(:,:,1)=I1;
% I(:,:,2)=I2;
% I(:,:,3)=I3;

Imgshow=uint8(Imgshow);
w=fspecial('gaussian',3,0.8);Imgshow=imfilter(Imgshow,w); 
imshow(Imgshow);
% figure(3),imshow(uint8(Imgshow(:,:,1)));
% figure(4),imshow(uint8(Imgshow(:,:,2)));
% figure(5),imshow(uint8(Imgshow(:,:,3)));










