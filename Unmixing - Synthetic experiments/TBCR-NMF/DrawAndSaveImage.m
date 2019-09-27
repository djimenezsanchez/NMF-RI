% 对比算法
BBmp = imread('bg2.bmp');

% maxV =6200;

% load('O:\开题报告和学位论文\论文数据\412老鼠数据\HALSL1\iter500sizehalf_lambda001.mat');
% CC = C.*(ones(Row*Col*scale2^2,1)*max(S,[],2)');
% % maxV = max(max(CC(:)))
% a = 2;
% figure,FluoImageShow3(reshape(CC(:,1),Row*scale2,Col*scale2),imresize(BBmp,scale2),0.01*max(CC(:,1))+10,maxV,a);
% figure,FluoImageShow3(reshape(CC(:,2),Row*scale2,Col*scale2),imresize(BBmp,scale2),0.01*max(CC(:,2))+10,maxV,a);
% figure,FluoImageShow3(reshape(CC(:,3),Row*scale2,Col*scale2),imresize(BBmp,scale2),0.01*max(CC(:,3))+10,maxV,a);
% % figure,FluoImageShow3(D4,BBmp,0.01*max(D4(:))+10,maxV,a);
% % figure,FluoImageShow3(D5,BBmp,0.01*max(D5(:))+15,maxV,a);
% 
% load('O:\开题报告和学位论文\论文数据\412老鼠数据\HALSsp\iter300sizehalfchannel5sp85.mat');
% CC = C.*(ones(Row*Col*scale2^2,1)*max(S,[],2)');
% figure,FluoImageShow3(reshape(CC(:,1),Row*scale2,Col*scale2),imresize(BBmp,scale2),0.01*max(CC(:,1))+10,maxV,a);
% figure,FluoImageShow3(reshape(CC(:,2),Row*scale2,Col*scale2),imresize(BBmp,scale2),0.01*max(CC(:,2))+10,maxV,a);
% figure,FluoImageShow3(reshape(CC(:,3),Row*scale2,Col*scale2),imresize(BBmp,scale2),0.01*max(CC(:,3))+10,maxV,a);
% % maxV = max(max(CC(:)))
% 
% load('O:\开题报告和学位论文\论文数据\412老鼠数据\SNMF\iter2000sizehalfspari85.mat');
% CC = C.*(ones(Row*Col*scale2^2,1)*max(S,[],2)');
% figure,FluoImageShow3(reshape(CC(:,1),Row*scale2,Col*scale2),imresize(BBmp,scale2),0.01*max(CC(:,1))+10,maxV,a);
% figure,FluoImageShow3(reshape(CC(:,2),Row*scale2,Col*scale2),imresize(BBmp,scale2),0.01*max(CC(:,2))+10,maxV,a);
% figure,FluoImageShow3(reshape(CC(:,3),Row*scale2,Col*scale2),imresize(BBmp,scale2),0.01*max(CC(:,3))+10,maxV,a);
% % maxV = max(max(CC(:)))

S = S./(max(S,[],2)*ones(1,4));
lw = 2;
figure,plot([525 542 579 624], S(1,:),'--bo','LineWidth',lw);
hold on;plot([525 542 579 624],S(2,:),'--cd','LineWidth',lw, 'Color',[0,0.5,0]);
hold on;plot([525 542 579 624],S(3,:),'--rs','LineWidth',lw);
axis([520 640 0 1.1]);
legend('AF488','AF555','Autofluorescence');
xlabel('Wavelength(nm)');
ylabel('Relative Intensity');

D = [D1(:) D2(:) D3(:) D4(:)];
OC = ((S*S')\(S*D'))';
OC(OC<0)=0;
temp = OC(:,1);
temp(temp < 0.02*max(temp))=0;
figure,plot(temp, 'b');xlabel('Pixel');ylabel('Abundance Intensity');axis([0 10^6 0 5800]);
temp = OC(:,2);
temp(temp < 0.02*max(temp))=0;
figure,plot(temp, 'g');xlabel('Pixel');ylabel('Abundance Intensity');axis([0 10^6 0 5800]);
temp = OC(:,3);
temp(temp < 0.02*max(temp))=0;
figure,plot(temp, 'r');xlabel('Pixel');ylabel('Abundance Intensity');axis([0 10^6 0 5800]);

scale2 = 1;


% load('O:\开题报告和学位论文\论文数据\1219日采集数据\SNMF\spa085.mat');
CC = C.*(ones(Row*Col*scale2^2,1)*max(S,[],2)');
maxV = 5800;
a = 2;
figure,FluoImageShow3(reshape(CC(:,1),Row*scale2,Col*scale2),imresize(BBmp,scale2),0.03*max(CC(:,1))+10,maxV,a);
figure,FluoImageShow3(reshape(CC(:,2),Row*scale2,Col*scale2),imresize(BBmp,scale2),0.03*max(CC(:,2))+10,maxV,a);
figure,FluoImageShow3(reshape(CC(:,3),Row*scale2,Col*scale2),imresize(BBmp,scale2),0.03*max(CC(:,3))+10,maxV,a);
% figure,FluoImageShow3(D4,BBmp,0.01*max(D4(:))+10,maxV,a);
% figure,FluoImageShow3(D5,BBmp,0.01*max(D5(:))+15,maxV,a);

load('O:\开题报告和学位论文\论文数据\1219日采集数据\SNMF\spa075.mat');
CC = C.*(ones(Row*Col*scale2^2,1)*max(S,[],2)');
maxV = max(max(CC(:)))
h = figure;FluoImageShow3(reshape(CC(:,1),Row*scale2,Col*scale2),imresize(BBmp,scale2),0.03*max(CC(:,1))+10,maxV,a);
% axis off;
% saveas(h,'O:\开题报告和学位论文\论文数据\1219日采集数据\SNMF\spa075-1','png');
h = figure;FluoImageShow3(reshape(CC(:,2),Row*scale2,Col*scale2),imresize(BBmp,scale2),0.03*max(CC(:,2))+10,maxV,a);
% axis off;
% saveas(h,'O:\开题报告和学位论文\论文数据\1219日采集数据\SNMF\spa075-2','png')
h = figure;FluoImageShow3(reshape(CC(:,3),Row*scale2,Col*scale2),imresize(BBmp,scale2),0.03*max(CC(:,3))+10,maxV,a);
% axis off;
% saveas(h,'O:\开题报告和学位论文\论文数据\1219日采集数据\SNMF\spa075-3','png')


load('O:\开题报告和学位论文\论文数据\1219日采集数据\HALSps\iter500sizehalfchannel5sp75.mat');
CC = C.*(ones(Row*Col*scale2^2,1)*max(S,[],2)');
maxV = max(max(CC(:)))
h = figure;FluoImageShow3(reshape(CC(:,1),Row*scale2,Col*scale2),imresize(BBmp,scale2),0.03*max(CC(:,1))+10,maxV,a);
% saveas(h,'O:\开题报告和学位论文\论文数据\1219日采集数据\SNMF\spa075-1','png');
h = figure;FluoImageShow3(reshape(CC(:,2),Row*scale2,Col*scale2),imresize(BBmp,scale2),0.03*max(CC(:,2))+10,maxV,a);
% saveas(h,'O:\开题报告和学位论文\论文数据\1219日采集数据\SNMF\spa075-2','png')
h = figure;FluoImageShow3(reshape(CC(:,3),Row*scale2,Col*scale2),imresize(BBmp,scale2),0.03*max(CC(:,3))+10,maxV,a);
% saveas(h,'O:\开题报告和学位论文\论文数据\1219日采集数据\SNMF\spa075-3','png')

load('O:\开题报告和学位论文\论文数据\1219日采集数据\HALSL1\iter500sizehalfchannel5lambda_001.mat');
CC = C.*(ones(Row*Col*scale2^2,1)*max(S,[],2)');
maxV = max(max(CC(:)))
h = figure;FluoImageShow3(reshape(CC(:,1),Row*scale2,Col*scale2),imresize(BBmp,scale2),0.03*max(CC(:,1))+10,maxV,a);
% saveas(h,'O:\开题报告和学位论文\论文数据\1219日采集数据\SNMF\spa075-1','png');
h = figure;FluoImageShow3(reshape(CC(:,2),Row*scale2,Col*scale2),imresize(BBmp,scale2),0.03*max(CC(:,2))+10,maxV,a);
% saveas(h,'O:\开题报告和学位论文\论文数据\1219日采集数据\SNMF\spa075-2','png')
h = figure;FluoImageShow3(reshape(CC(:,3),Row*scale2,Col*scale2),imresize(BBmp,scale2),0.03*max(CC(:,3))+10,maxV,a);

close all;