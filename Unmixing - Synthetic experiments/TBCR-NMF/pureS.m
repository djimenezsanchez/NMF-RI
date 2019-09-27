close all;clear all;
% 标定试验
BBmp = imread('96孔板标定试验\474-579.bmp');
data = load('96孔板标定试验\474-579.txt');
data = medfilt2(data);
data= FluoNoiseReduction(data);
[m,n] = size(data);

% 去除下面一行
% data(400:800,:)=0;

maxV = max(data(:));
% figure,FluoImageShow3(data,BBmp, 0.01*maxV,maxV,2);
% figure, imagesc(data);

% 去除一些不是目标的区域
mask = false(m,n);
mask(data > 0.02*maxV) = true;
mask2 = imopen(mask,strel('disk',2));
mask2(1:200,:) = false;
mask2(:,940:end) = false;
figure,imshow(mask2);

hole488 = false(m,n);
hole488(1:400,380:580) = mask2(1:400,380:580);%figure,imshow(m13);
hole488 = imerode(hole488, strel('disk',2));
figure,imshow(hole488);

hole555 = false(m,n);
hole555(560:780,380:550) = mask2(560:780,380:550);
hole555 = imerode(hole555,strel('disk',2));
figure,imshow(hole555);

data1 = load('96孔板标定试验\474-525.txt');data1 = medfilt2(data1);data1= FluoNoiseReduction(data1);
data2 = load('96孔板标定试验\474-542.txt');data2 = medfilt2(data2);data2= FluoNoiseReduction(data2);
data3 = load('96孔板标定试验\500-579.txt');data3 = medfilt2(data3);data3= FluoNoiseReduction(data3);
data4 = load('96孔板标定试验\500-624.txt');data4 = medfilt2(data4);data4= FluoNoiseReduction(data4);
s488 = [mean(data1(hole488)) mean(data2(hole488)) mean(data3(hole488)) mean(data4(hole488))];
s555 = [mean(data1(hole555)) mean(data2(hole555)) mean(data3(hole555)) mean(data4(hole555)) ];

bg1 = load('老鼠背景\474-525.txt');bg1 = medfilt2(bg1);bg1 = FluoNoiseReduction(bg1);
bg2 = load('老鼠背景\474-542.txt');bg2 = medfilt2(bg2);bg2 = FluoNoiseReduction(bg2);
bg3 = load('老鼠背景\500-579.txt');bg3 = medfilt2(bg3);bg3 = FluoNoiseReduction(bg3);
bg4 = load('老鼠背景\500-624.txt');bg4 = medfilt2(bg4);bg4 = FluoNoiseReduction(bg4);

bgmask = false(m,n);
bgmask(bg1 > 0.02*max(bg1(:))) = true;
bgmask = imerode(bgmask,strel('disk',2));
figure,imshow(bgmask);
sAF = [mean(bg1(bgmask)) mean(bg2(bgmask)) mean(bg3(bgmask)) mean(bg4(bgmask))];

s488 = s488/max(s488(:));
s555 = s555/max(s555(:));
sAF = sAF/max(sAF(:));
lw = 2;
figure,plot([525 542 579 624], s488,'--bo','LineWidth',lw);
hold on;plot([525 542 579 624],s555,'--cd','LineWidth',lw, 'Color',[0,0.5,0]);
hold on;plot([525 542 579 624],sAF,'--rs','LineWidth',lw);
axis([520 640 0 1.1]);
legend('AF488','AF555','Autofluorescence');
xlabel('Wavelength(nm)');
ylabel('Relative Intensity');
% data1 = load('96孔板标定试验\474-525.txt');data1 = medfilt2(data1);data1= FluoNoiseReduction(data1);
% data2 = load('96孔板标定试验\474-542.txt');data2 = medfilt2(data2);data2= FluoNoiseReduction(data2);
% data3 = load('96孔板标定试验\474-579.txt');data3 = medfilt2(data3);data3= FluoNoiseReduction(data3);
% data4 = load('96孔板标定试验\500-624.txt');data4 = medfilt2(data4);data4= FluoNoiseReduction(data4);
% data5 = load('96孔板标定试验\500-716.txt');data5 = medfilt2(data5);data5= FluoNoiseReduction(data5);
% s488 = [mean(data1(hole488)) mean(data2(hole488)) mean(data3(hole488)) mean(data4(hole488)) mean(data5(hole488)) ];
% s555 = [mean(data1(hole555)) mean(data2(hole555)) mean(data3(hole555)) mean(data4(hole555)) mean(data5(hole555)) ];
% 
% bg1 = load('老鼠背景\474-525.txt');bg1 = medfilt2(bg1);bg1 = FluoNoiseReduction(bg1);
% bg2 = load('老鼠背景\474-542.txt');bg2 = medfilt2(bg2);bg2 = FluoNoiseReduction(bg2);
% bg3 = load('老鼠背景\474-579.txt');bg3 = medfilt2(bg3);bg3 = FluoNoiseReduction(bg3);
% bg4 = load('老鼠背景\500-624.txt');bg4 = medfilt2(bg4);bg4 = FluoNoiseReduction(bg4);
% bg5 = load('老鼠背景\500-716.txt');bg5 = medfilt2(bg5);bg5 = FluoNoiseReduction(bg5);
% bgmask = false(m,n);
% bgmask(bg1 > 0.02*max(bg1(:))) = true;
% bgmask = imerode(bgmask,strel('disk',2));
% figure,imshow(bgmask);
% sAF = [mean(bg1(bgmask)) mean(bg2(bgmask)) mean(bg3(bgmask)) mean(bg4(bgmask)) mean(bg5(bgmask)) ];
% 
% s488 = s488/max(s488(:));
% s555 = s555/max(s555(:));
% sAF = sAF/max(sAF(:));
% figure,plot(1:5, s488,'b','LineWidth',3);
% hold on;plot(1:5,s555,'g','LineWidth',3);
% hold on;plot(1:5,sAF,'r','LineWidth',3);
% 
% load('O:\开题报告和学位论文\论文数据\1219日采集数据\HALSps\iter500sizehalfchannel5sp85.mat', 'S');
% load('O:\开题报告和学位论文\论文数据\1219日采集数据\SNMF\spa075.mat', 'S');
% load('O:\开题报告和学位论文\论文数据\1219日采集数据\HALSL1\iter500sizehalfchannel5lambda_001.mat', 'S');
% os = S;
% os(1,:) = os(1,:)/max(os(1,:));
% os(2,:) = os(2,:)/max(os(2,:));
% os(3,:) = os(3,:)/max(os(3,:));
% figure,plot(os','LineWidth',3);
% 
% SADALL = SAD_3D([s488;s555;sAF], os)
% mean(SADALL)
