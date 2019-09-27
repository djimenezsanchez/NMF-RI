% Normlized Cut
close all;clear all;

BBmp = imread('bg.bmp');
% D1 = load('threefluo\474-525.txt');
% D2 = load('threefluo\474-542.txt');
% D3 = load('threefluo\474-579.txt');
% D4 = load('threefluo\474-624.txt');
% D5 = load('threefluo\474-716.txt');
% D1 = load('新数据\474-525.txt');
% D2 = load('新数据\474-542.txt');
% D3 = load('新数据\474-579.txt');
% D4 = load('新数据\474-624.txt');
% D5 = load('新数据\474-716.txt');
% D1 = load('新数据2\474-525.txt');
% D2 = load('新数据2\474-542.txt');
% D3 = load('新数据2\474-579.txt');
% D4 = load('新数据2\474-624.txt');
% D5 = load('新数据2\474-716.txt');
% D1 = load('474-525(3).txt');
% D2 = load('474-542(3).txt');
% D3 = load('474-579(3).txt');
% D4 = load('500-624(3).txt');
% D5 = load('500-716(3).txt');
% D1 = load('onefluo\474-525.txt');
% D2 = load('onefluo\474-542.txt');
% D3 = load('onefluo\474-579.txt');
% D4 = load('onefluo\474-624.txt');
% D5 = load('onefluo\474-716.txt');
D1 = load('fourfluo\474-525.txt');
D2 = load('fourfluo\474-542.txt');
D3 = load('fourfluo\500-579.txt');
D4 = load('fourfluo\500-624.txt');
D5 = load('fourfluo\500-716.txt');
% D1 = load('zy\525.txt');D1 = D1(1:1000,1:1000);
% D2 = load('zy\542.txt');D2 = D2(1:1000,1:1000);
% D3 = load('zy\579.txt');D3 = D3(1:1000,1:1000);
% D4 = load('zy\624.txt');D4 = D4(1:1000,1:1000);
% D5 = load('zy\716.txt');D5 = D5(1:1000,1:1000);

% D1 = load('裸鼠\474-525.txt');
% D2 = load('裸鼠\474-542.txt');
% D3 = load('裸鼠\474-579.txt');
% D4 = load('裸鼠\474-624.txt');
% D5 = load('裸鼠\474-716.txt');

figure,imagesc([D1 D2 D3 D4 D5]);

% 滤波以及去除基底信号
% D1 = medfilt2(D1);
% D2 = medfilt2(D2);
% D3 = medfilt2(D3);
% D4 = medfilt2(D4);
% D5 = medfilt2(D5);
D1 = FluoNoiseReduction(D1);
D2 = FluoNoiseReduction(D2);
D3 = FluoNoiseReduction(D3);
D4 = FluoNoiseReduction(D4);
D5 = FluoNoiseReduction(D5);


maxV = max(max([D1 D2 D3 D4 D5]));
a = 2;
figure,FluoImageShow3(D1,BBmp,0.01*max(D1(:))+10,maxV,a);
figure,FluoImageShow3(D2,BBmp,0.01*max(D2(:))+10,maxV,a);
figure,FluoImageShow3(D3,BBmp,0.01*max(D3(:))+10,maxV,a);
figure,FluoImageShow3(D4,BBmp,0.01*max(D4(:))+10,maxV,a);
figure,FluoImageShow3(D5,BBmp,0.01*max(D5(:))+15,maxV,a);

% temp = zeros(1000,1000);
% if size(D1,1) == 992
%     temp(1:992,1:992) = D1;D1 = temp;
%     temp(1:992,1:992) = D2;D2 = temp;
%     temp(1:992,1:992) = D3;D3 = temp;
%     temp(1:992,1:992) = D4;D4 = temp;
%     temp(1:992,1:992) = D5;D5 = temp;
% end

% tar1 = imread('仿真数据1\target1.jpg');
% tar2 = imread('仿真数据1\target2.jpg');
% auto = imread('仿真数据1\auto.jpg');
% % tar1 = imread('仿真数据2\AF594V1.bmp');
% % tar2 = imread('仿真数据2\AF488V1.bmp');
% % auto = imread('仿真数据2\background1.bmp');
% tar1 = double(tar1)/double(max(tar1(:)));
% tar2 = double(tar2)/double(max(tar2(:)));
% auto = double(auto)/double(max(auto(:)));
% 
% st1 = [0.2  0.3     0.8     1        0.4];
% st2 = [1	 0.8     0.6     0.1     0];
% sa = [0.5   1        0.8     0.6     0.6];
% figure,plot(1:5,st1,1:5,st2,1:5,sa);
% 
% D1 = tar1*st1(1) + tar2*st2(1) + auto*sa(1);
% D2 = tar1*st1(2) + tar2*st2(2) + auto*sa(2);
% D3 = tar1*st1(3) + tar2*st2(3) + auto*sa(3);
% D4 = tar1*st1(4) + tar2*st2(4) + auto*sa(4);
% D5 = tar1*st1(5) + tar2*st2(5) + auto*sa(5);


[Row,Col] = size(D1);
% data downsampling for the subsequent Normalized Cut
if Row > 500
    scale = 0.1;
else
    scale = 0.25;
end
d1 = imresize(D1, scale);
d2 = imresize(D2, scale);
d3 = imresize(D3, scale);
d4 = imresize(D4, scale);
d5 = imresize(D5, scale);
[row, col] = size(d1);

% compute mask
mask = false(row,col);
mask(d1 > 0.02*max(d1(:))) = true;
mask = imerode(mask,strel('disk',1));
m = sum(mask(:));
X = zeros(m,5);
X(:,1) = d1(mask);
X(:,2) = d2(mask);
X(:,3) = d3(mask);
X(:,4) = d4(mask);
X(:,5) = d5(mask);

%%
% 坐标信息
% [CoordX CoordY] = find(mask == true);
% CoordX = zeros(m,1);
% CoordY = zeros(m,1);
% mask2 = false(row,col);
% count = 1;
% for i = 1:row
%     for j = 1:col
%         if mask(i,j) == true
%             CoordX(count) = i;
%             CoordY(count) = j;
%             count = count + 1;
%         end
%      end
% end
% 
% for i = 1:m
%         mask2(CoordX(i),CoordY(i)) = true;
% end
% figure,imshow(mask2);
%%
orignalmask = mask;                    % 保存原始的整个老鼠的掩膜
[y,v,nv,evec,eval] = Ncut(X);nv
res = zeros(row,col);
res(mask) = v;
% figure,imagesc(res);

binpar = false(row,col);
binpar(res>0.0) = true;
% figure,imshow(binpar);


seg = zeros(row,col);
reg = 1;
seg(mask) = reg;
reg = reg + 1;
if sum(binpar(:)) > sum(mask(:))/2
    seg(mask & ~binpar) = reg;
    mask = mask & binpar;
else
    seg(mask & binpar) = reg;
    mask = mask & ~binpar;
end

NcutValue = 0.2;
while eval(2,2) < NcutValue
%     if sum(binpar(:)) > sum(mask(:))/2
%         mask = mask & binpar;
%     else
%         mask = mask & ~binpar;
%     end
    m = sum(mask(:)); 
    X = zeros(m,5);
    X(:,1) = d1(mask);
    X(:,2) = d2(mask);
    X(:,3) = d3(mask);
    X(:,4) = d4(mask);
    X(:,5) = d5(mask);
    [y,v,nv,evec,eval] = Ncut(X); nv
    if eval(2,2) < NcutValue
        res = zeros(row,col);
        res(mask) = v;
        figure,imagesc(res);
        binpar = false(row,col);
        binpar(res>0.0) = true;
%         figure,imshow(binpar);

        reg = reg + 1;
        if sum(binpar(:)) > sum(mask(:))/2
            seg(mask & ~binpar) = reg;
            mask = mask & binpar;
        else
            seg(mask & binpar) = reg;
            mask = mask & ~binpar;
        end
    end
end

figure,imagesc(seg);
% figure,imshow(mask);

%% 端元提取,其这mask此时此刻已经是纯自发荧光的mask了,targetmask为目标的mask.然后将目标区域减去背景的平均值然后再运用ATGP进行端元提取
targetmask = xor(orignalmask , mask);
figure,imshow(targetmask);
TX = zeros(sum(targetmask(:)), 5);
TX(:,1) = d1(targetmask);
TX(:,2) = d2(targetmask);
TX(:,3) = d3(targetmask);
TX(:,4) = d4(targetmask);
TX(:,5) = d5(targetmask);
AverageAuto = [mean(d1(mask)) mean(d2(mask)) mean(d3(mask)) mean(d4(mask)) mean(d5(mask))];
TX = TX - repmat(AverageAuto,sum(targetmask(:)),1);
TX(TX < 0) = 0;

% [E,index]=ATGP(TX,2);

[E,C] = EIA_ATGP(TX',2);
E = E./ repmat(max(E),5,1);
figure,plot(1:5,E(:,1),1:5,E(:,2),1:5,AverageAuto/max(AverageAuto));

% figure,plot(1:5,st1,1:5,st2,1:5,sa);
% [E, indice, Rp] = VCA(TX','Endmembers',3);
% E(E<0)=0;
% E = E./ repmat(max(E),5,1);
% figure,plot(1:5,E(:,1),1:5,E(:,2),1:5,AverageAuto/max(AverageAuto));

%% 接下来就是对原始数据进行分解,这个时候的方法很多
scale2 = 0.5;eps = 1e-9;
V = zeros(Row * Col*scale2*scale2, 5);
SD1 = imresize(D1, scale2);SD1(SD1<eps) = eps;
SD2 = imresize(D2, scale2);SD2(SD2<eps) = eps;
SD3 = imresize(D3, scale2);SD3(SD3<eps) = eps;
SD4 = imresize(D4, scale2);SD4(SD4<eps) = eps;
SD5 = imresize(D5, scale2);SD5(SD5<eps) = eps;
V(:,1) = SD1(:);
V(:,2) = SD2(:);
V(:,3) = SD3(:);
V(:,4) = SD4(:);
V(:,5) = SD5(:);
phi = 0.80; targetNum = 2;lambda = 0.01;
A0 = rand(Row*Col, targetNum + 1);
S0 = rand(targetNum + 1,5);
S0(1:2,:) = E';
S0(targetNum + 1, :) = AverageAuto/max(AverageAuto);
[C0,S0]=Initial_Pure(V,S0, 100);

phi = 0.85;
[C,S,time,cost_record,lambda_sq]= mynmf_ghals_L12_Ultimate(V,C0,S0,phi);
% % [C,S]= myHALS12(V,C0,S0,phi);
% % [C,S,time,cost_record]= mynmf_ghals(V,C0,S0,lambda);
% % [C,S,t,objhistory]=myNMFsv(V,phi,C0,S0);%,tol)
% 
% % [C,S]=myNMFMU(V, S0,2000);
% figure,imagesc(reshape(C(:,1),Row*scale2,Col*scale2));
% figure,imagesc(reshape(C(:,2),Row*scale2,Col*scale2));
% figure,imagesc(reshape(C(:,3),Row*scale2,Col*scale2));
% figure,plot(S');
% figure,imagesc(reshape(C0(:,1),Row*scale2,Col*scale2));
% figure,imagesc(reshape(C0(:,2),Row*scale2,Col*scale2));
% figure,imagesc(reshape(C0(:,3),Row*scale2,Col*scale2));


