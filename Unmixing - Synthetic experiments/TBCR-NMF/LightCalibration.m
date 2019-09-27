close all;clear all;
% 标定试验
BBmp = imread('96孔板标定试验\474-579.bmp');
data = load('96孔板标定试验\474-525.txt');
data = medfilt2(data);
data= FluoNoiseReduction(data);
[m,n] = size(data);

% 去除下面一行
data(400:800,:)=0;

maxV = max(data(:));
figure,FluoImageShow3(data,BBmp, 0.01*maxV,maxV,2);
figure, imagesc(data);

% 去除一些不是目标的区域
mask = false(m,n);
mask(data > 0.02*maxV) = true;
mask2 = imopen(mask,strel('disk',2));
mask2(1:200,:) = false;
mask2(:,940:end) = false;
% figure,imshow(mask2);

m11 = false(m,n);
m11(1:400,1:180) = mask2(1:400,1:180); %figure,imshow(m11);
m12 = false(m,n);
m12(1:400,180:380) = mask2(1:400,180:380); %figure,imshow(m12);
m13 = false(m,n);
m13(1:400,380:580) = mask2(1:400,380:580);%figure,imshow(m13);
m14 = false(m,n);
m14(1:400,580:780) = mask2(1:400,580:780); %figure,imshow(m14);
m15 = false(m,n);
m15(1:400,780:end) = mask2(1:400,780:end);%figure,imshow(m15);

m21 = false(m,n);
m21(400:800,1:180) = mask2(400:800,1:180);% figure,imshow(m21);
m22 = false(m,n);
m22(400:800,180:380) = mask2(400:800,180:380); %figure,imshow(m22);
m23 = false(m,n);
m23(400:800,380:580) = mask2(400:800,380:580); %figure,imshow(m23);
m24 = false(m,n);
m24(400:800,580:780) = mask2(400:800,580:780); %figure,imshow(m24);
m25 = false(m,n);
m25(400:800,780:end) = mask2(400:800,780:end);% figure,imshow(m25);

lightdis = load('96孔板标定试验\lightdis3.txt');
lightdis = imfilter(lightdis, fspecial('gaussian',[3 3]));
% lightdis = imfilter(lightdis, fspecial('unsharp'));
lightdis = medfilt2(lightdis);
lightdis = imfilter(lightdis, fspecial('average',3));
lightdis = lightdis - 303;
lightdis(lightdis < 0) = 0;
% remask = false(m,n);
% remask(lightdis > 1000) = true;
% figure,imshow(remask);
% lightdis = roifill(lightdis, remask);
% figure,imagesc(lightdis);axis image;axis off;



targetdis = ones(m,n)*1000;
targetdis(mask2) = lightdis(mask2);
% figure,imagesc(targetdis);
Ntargetdis = targetdis/median(targetdis(mask));
figure,imagesc(Ntargetdis);

cali_data = data./Ntargetdis;
figure,imagesc(cali_data);


oa= [mean(data(m11)) ...
      mean(data(m12)) ...
      mean(data(m13)) ...
      mean(data(m14)) ...
      mean(data(m15))]
% ob = [mean(data(m21)) ...
%        mean(data(m22)) ...
%        mean(data(m23)) ...
%        mean(data(m24)) ...
%        mean(data(m25))]

a= [mean(cali_data(m11)) ...
      mean(cali_data(m12)) ...
      mean(cali_data(m13)) ...
      mean(cali_data(m14)) ...
      mean(cali_data(m15))]
% b = [mean(cali_data(m21)) ...
%        mean(cali_data(m22)) ...
%        mean(cali_data(m23)) ...
%        mean(cali_data(m24)) ...
%        mean(cali_data(m25))]
% 去除坏点
% relightdis = roifill(lightdis,remask);
% figure,imagesc(relightdis);
maxV = max(cali_data(:));
figure,FluoImageShow3(cali_data,BBmp, 0.01*maxV,maxV,2);







