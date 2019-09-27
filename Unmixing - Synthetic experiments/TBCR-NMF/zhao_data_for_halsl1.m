close all;
% 赵勇数据计算
% D1 = load('有毛鼠数据\525.txt');%D1 = D1(1:1000,1:1000);
% D2 = load('有毛鼠数据\542.txt');%D2 = D2(1:1000,1:1000);
% D3 = load('有毛鼠数据\579.txt');%D3 = D3(1:1000,1:1000);
% D4 = load('有毛鼠数据\624.txt');%D4 = D4(1:1000,1:1000);
% D5 = load('有毛鼠数据\716.txt');%D5 = D5(1:1000,1:1000);
D1 = load('无毛鼠数据\5252.txt');%D1 = D1(1:1000,1:1000);
D2 = load('无毛鼠数据\5422.txt');%D2 = D2(1:1000,1:1000);
D3 = load('无毛鼠数据\5792.txt');%D3 = D3(1:1000,1:1000);
D4 = load('无毛鼠数据\6242.txt');%D4 = D4(1:1000,1:1000);
D5 = load('无毛鼠数据\7162.txt');%D5 = D5(1:1000,1:1000);
figure,imagesc([D1 D2 D3 D4 D5]);

% 滤波以及去除基底信号
D1 = medfilt2(D1);
D2 = medfilt2(D2);
D3 = medfilt2(D3);
D4 = medfilt2(D4);
D5 = medfilt2(D5);
D1 = FluoNoiseReduction(D1);
D2 = FluoNoiseReduction(D2);
D3 = FluoNoiseReduction(D3);
D4 = FluoNoiseReduction(D4);
D5 = FluoNoiseReduction(D5);

D1 = imresize(D1,0.5);
D2 = imresize(D2,0.5);
D3 = imresize(D3,0.5);
D4 = imresize(D4,0.5);
D5 = imresize(D5,0.5);

[Row Col] = size(D1);
V = zeros(Row * Col, 5);
V(:,1) = D1(:);
V(:,2) = D2(:);
V(:,3) = D3(:);
V(:,4) = D4(:);
V(:,5) = D5(:);

S0 = rand(3, 5);
[C0,S0]=Initial_Pure(V,S0,200);
% C0 = rand(Row*Col, 3);

lambda = 0.01;
 [C,S,time,cost_record]= mynmf_ghals(V,C0,S0,lambda);
 
% figure,imagesc(reshape(C(:,1),Row,Col));
% figure,imagesc(reshape(C(:,2),Row,Col));
% figure,imagesc(reshape(C(:,3),Row,Col));


maxV = max(C(:));
a = 20;
th = 0.01;
close all;
BBmp = imread('无毛鼠数据\BG.bmp');
figure,FluoImageShow3(reshape(C(:,1),Row,Col),imresize(BBmp,0.5),th*max(C(:,1)),maxV,a);
figure,FluoImageShow3(reshape(C(:,2),Row,Col),imresize(BBmp,0.5),th*max(C(:,2)),maxV,a);
BBmp = imread('无毛鼠数据\BG2.bmp');
figure,FluoImageShow3(reshape(C(:,3),Row,Col),imresize(BBmp,0.5),th*max(C(:,3)),maxV,a);
