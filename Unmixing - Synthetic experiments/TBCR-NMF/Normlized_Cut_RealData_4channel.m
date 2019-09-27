% Normlized Cut
% close all;clear all;


D1 = load('imagedata\474-525.txt');
D2 = load('imagedata\474-542.txt');
D3 = load('imagedata\500-579.txt');
D4 = load('imagedata\500-624.txt');

% D1 = load('threefluo\474-525.txt');
% D2 = load('threefluo\474-542.txt');
% D3 = load('threefluo\500-579.txt');
% D4 = load('threefluo\500-624.txt');

%% 
% D1 = load('474-525.txt');
% D2 = load('474-542.txt');
% D3 = load('474-579.txt');
% D4 = load('474-624.txt');
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
% D1 = load('fourfluo\474-525.txt');
% D2 = load('fourfluo\474-542.txt');
% D3 = load('fourfluo\500-579.txt');
% D4 = load('fourfluo\500-624.txt');
% D5 = load('fourfluo\474-716.txt');
% D1 = load('zy\525.txt');D1 = D1(1:1000,1:1000);
% D2 = load('zy\542.txt');D2 = D2(1:1000,1:1000);
% D3 = load('zy\579.txt');D3 = D3(1:1000,1:1000);
% D4 = load('zy\624.txt');D4 = D4(1:1000,1:1000);
% D5 = load('zy\716.txt');D5 = D5(1:1000,1:1000);

% D1 = load('mouse\474-525.txt');
% D2 = load('mouse\474-542.txt');
% D3 = load('mouse\474-579.txt');
% D4 = load('mouse\474-624.txt');
% D5 = load('mouse\474-716.txt');
%%



figure,imagesc([D1 D2 D3 D4]);

% Filter out the noise and baseline signal
D1 = medfilt2(D1);
D2 = medfilt2(D2);
D3 = medfilt2(D3);
D4 = medfilt2(D4);

D1 = FluoNoiseReduction(D1);
D2 = FluoNoiseReduction(D2);
D3 = FluoNoiseReduction(D3);
D4 = FluoNoiseReduction(D4);


BBmp = imread('imagedata\bg.bmp');
% BBmp = imread('newdata3\bg.bmp');
% BBmp = imread('fourfluo\bg.bmp');
maxV = max(max([D1 D2 D3 D4]));
a = 2;
figure,FluoImageShow3(D1,BBmp,0.01*max(D1(:))+10,maxV,a);
figure,FluoImageShow3(D2,BBmp,0.01*max(D2(:))+10,maxV,a);
figure,FluoImageShow3(D3,BBmp,0.01*max(D3(:))+10,maxV,a);
figure,FluoImageShow3(D4,BBmp,0.01*max(D4(:))+10,maxV,a);

% D1 = imresize(D1, 0.5);
% D2 = imresize(D2, 0.5);
% D3 = imresize(D3, 0.5);
% D4 = imresize(D4, 0.5);

% temp = zeros(1000,1000);
% if size(D1,1) == 992
%     temp(1:992,1:992) = D1;D1 = temp;
%     temp(1:992,1:992) = D2;D2 = temp;
%     temp(1:992,1:992) = D3;D3 = temp;
%     temp(1:992,1:992) = D4;D4 = temp;
%     temp(1:992,1:992) = D5;D5 = temp;
% end

% tar1 = imread('simulation1\target1.jpg');
% tar2 = imread('simulation1\target2.jpg');
% auto = imread('simulation1\auto.jpg');
% % tar1 = imread('simulation2\AF594V1.bmp');
% % tar2 = imread('simulation2\AF488V1.bmp');
% % auto = imread('simulation2\background1.bmp');
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
% Data down-sampling for subsequent Normalized Cut
if Row > 500
    scale = 0.1;
else
    scale = 0.25;
end

d1 = imresize(D1, scale);
d2 = imresize(D2, scale);
d3 = imresize(D3, scale);
d4 = imresize(D4, scale);
[row, col] = size(d1);

% Mask computing
mask = false(row,col);
mask(d1 > 0.01*max(d1(:))) = true;
mask = imerode(mask,strel('disk',1));
m = sum(mask(:));
X = zeros(m,4);
X(:,1) = d1(mask);
X(:,2) = d2(mask);
X(:,3) = d3(mask);
X(:,4) = d4(mask);

%%
tic;
orignalmask = mask;                    % the mask of whole mouse
[y,v,nv,evec,eval] = Ncut(X);%nv
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

NcutValue = 0.4;
while eval(2,2) < NcutValue
%     if sum(binpar(:)) > sum(mask(:))/2
%         mask = mask & binpar;
%     else
%         mask = mask & ~binpar;
%     end
    eval(2,2)
    m = sum(mask(:)); 
    X = zeros(m,4);
    X(:,1) = d1(mask);
    X(:,2) = d2(mask);
    X(:,3) = d3(mask);
    X(:,4) = d4(mask);

    [y,v,nv,evec,eval] = Ncut(X);%nv
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

%% Endmember extraction using ATGP after the average background autofluorescence are removed from the target regions
%% the mask is the mask from the background autofluorescence¡À,targetmask is the mask for the target regions.
targetmask = xor(orignalmask , mask);
figure,imshow(targetmask);
TX = zeros(sum(targetmask(:)), 4);
TX(:,1) = d1(targetmask);
TX(:,2) = d2(targetmask);
TX(:,3) = d3(targetmask);
TX(:,4) = d4(targetmask);
AverageAuto = [mean(d1(mask)) mean(d2(mask)) mean(d3(mask)) mean(d4(mask))];
TX = TX - repmat(AverageAuto,sum(targetmask(:)),1);
TX(TX < 0) = 0;

[E,C] = EIA_ATGP(TX',2);
E = E./ repmat(max(E),4,1);
figure,plot(1:4,E(:,1),1:4,E(:,2),1:4,AverageAuto/max(AverageAuto));
toc;
% 
% [E,C] = EIA_ATGP(TX',3);
% E = E./ repmat(max(E),4,1);
% figure,plot(1:4,E(:,1),1:4,E(:,2),1:4,E(:,3),1:4,AverageAuto/max(AverageAuto));
% figure,plot(1:5,st1,1:5,st2,1:5,sa);
% [E, indice, Rp] = VCA(TX','Endmembers',2);
% % E(E<0)=0;
% E = E./ repmat(max(E),4,1);
% figure,plot(1:4,E(:,1),1:4,E(:,2),1:4,AverageAuto/max(AverageAuto));

%% unmixing the original mixed data
bin = 1;
scale2 = 1/bin;
V = zeros(Row * Col / (bin*bin), 4);
% SD1 = imresize(D1, scale2);SD1(SD1<eps) = eps;
% SD1 = DoBinning2(D1);
% SD2 = DoBinning2(D2);
% SD3 = DoBinning2(D3);
% SD4 = DoBinning2(D4);
% V(:,1) = SD1(:);
% V(:,2) = SD2(:);
% V(:,3) = SD3(:);
% V(:,4) = SD4(:);
V(:,1) = D1(:);
V(:,2) = D2(:);
V(:,3) = D3(:);
V(:,4) = D4(:);
phi = 0.85; targetNum = 2;
A0 = rand(Row*Col, targetNum + 1);
S0 = rand(targetNum + 1,4);
S0(1:targetNum,:) = E'; % initializing the target fluorescences using the above-extracted target endmember
S0(targetNum + 1, :) = AverageAuto/max(AverageAuto); % the background autofluorescece
%  [C,S]=myNMFMU(V, S0,2000);
[C0,S0]=Initial_Pure(V,S0,100);

tic;
phi = 0.80; 
[C,S,time,cost_record,lambda_sq]= mynmf_ghals_L12_Ultimate(V,C0,S0,phi);
toc;
% save('D:\MATLABProject\mouse3\HALS_L12_phi080.mat', 'C','S');
% phi = 0.75;
% [C,S,time,cost_record,lambda_sq]= mynmf_ghals_L12_Ultimate(V,C0,S0,phi);
% save('D:\MATLABProject\mouse3\HALS_L12_phi075.mat', 'C','S');
% 
% [C,S,e] = LagrangianNMU(V,3,1000,C0,S0);                  %NMU
% save('D:\MATLABProject\mouse3\NMU1000.mat', 'C','S');
% 
% maskL = false(size(SD1));
% mask(SD1 > 0.01*max(SD1(:))) = true;
% L = sum(mask(:));
% [C, S, errorHistory] = sparseNMF_W(V, 3, L, 300, 10, 1); %NMFL0
% save('D:\MATLABProject\mouse3\NMFL0_iter300_10.mat', 'C','S');
% [C, S, errorHistory] = sparseNMF_W(V, 3, L, 200, 30, 1); %NMFL0
% save('D:\MATLABProject\mouse3\NMFL0_iter200_30.mat', 'C','S');
% 
% tic;
% lambda = 20;
% [C,S,time,cost_record]= mynmf_ghals(V,C0,S0,lambda);
% toc;
% save('D:\MATLABProject\mouse3\HALS_L1_lambda002.mat', 'C','S');
% lambda = 0.01;
% [C,S,time,cost_record]= mynmf_ghals(V,C0,S0,lambda);
% save('D:\MATLABProject\mouse3\HALS_L1_lambda001.mat', 'C','S');
% lambda = 0.005;
% [C,S,time,cost_record]= mynmf_ghals(V,C0,S0,lambda);
% save('D:\MATLABProject\mouse3\HALS_L1_lambda0005.mat', 'C','S');
% 
% 
% phi = 0.75;
% [C,S,time,cost_record]= mynmf_ghals_L1(V,C0,S0,phi);
% save('D:\MATLABProject\mouse3\HALS_L1_phi075.mat', 'C','S');
% phi = 0.80;
% [C,S,time,cost_record]= mynmf_ghals_L1(V,C0,S0,phi);
% save('D:\MATLABProject\mouse3\HALS_L1_phi080.mat', 'C','S');
% 
% 
% tic;
% phi = 0.70;
% [C,S,t,objhistory]=myNMFsv(V,phi,C0,S0);
% % save('D:\MATLABProject\mouse3\NMFsv_phi075.mat', 'C','S');
% toc;
% phi = 0.80;
% [C,S,t,objhistory]=myNMFsv(V,phi,C0,S0);
% save('D:\MATLABProject\mouse3\NMFsv_phi080.mat', 'C','S');
% phi = 0.85;
% [C,S,t,objhistory]=myNMFsv(V,phi,C0,S0);
% save('D:\MATLABProject\mouse3\NMFsv_phi085.mat', 'C','S');

%% displayAll
    S = S./(max(S,[],2)*ones(1,4));
    lw = 2;
    figure,plot([525 542 579 624], S(1,:),'--bo','LineWidth',lw);
    hold on;plot([525 542 579 624],S(2,:),'--cd','LineWidth',lw, 'Color',[0,0.5,0]);
    hold on;plot([525 542 579 624],S(3,:),'--rs','LineWidth',lw);
    axis([520 640 0 1.1]);
    legend('AF488','AF555','Autofluorescence');
    xlabel('Wavelength(nm)');
    ylabel('Relative Intensity');
% 
%     D = [D1(:) D2(:) D3(:) D4(:)];
%     OC = ((S*S')\(S*D'))';
%     OC(OC<0)=0;
%     temp = OC(:,1);
%     temp(temp < 0.02*max(temp))=0;
%     figure,plot(temp, 'b');xlabel('Pixel');ylabel('Abundance Intensity');axis([0 10^6 0 5800]);
%     temp = OC(:,2);
%     temp(temp < 0.02*max(temp))=0;
%     figure,plot(temp, 'g');xlabel('Pixel');ylabel('Abundance Intensity');axis([0 10^6 0 5800]);
%     temp = OC(:,3);
%     temp(temp < 0.02*max(temp))=0;
%     figure,plot(temp, 'r');xlabel('Pixel');ylabel('Abundance Intensity');axis([0 10^6 0 5800]);
% 
% 
%     % load('O:\MasterProject\1219\SNMF\spa085.mat');
    CC = C.*(ones(Row*Col*scale2^2,1)*max(S,[],2)');
    maxV = max(CC(:));
    a = 2;
    figure,FluoImageShow3(reshape(CC(:,1),Row*scale2,Col*scale2),imresize(BBmp,scale2),0.03*max(CC(:,1)),maxV,a);
    figure,FluoImageShow3(reshape(CC(:,2),Row*scale2,Col*scale2),imresize(BBmp,scale2),0.03*max(CC(:,2)),maxV,a);
    figure,FluoImageShow3(reshape(CC(:,3),Row*scale2,Col*scale2),imresize(BBmp,scale2),0.03*max(CC(:,3)),maxV,a);
% SADALL = SAD_3D([s488;s555;sAF], S)
% close all;
% [C,S]=myNMFMU(V,C,S,2000);
