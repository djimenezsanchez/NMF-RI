for spa=65:5:85
close all;clearvars -except spa;
x1 = 480;x2=650;step=20;
[s488 s555 sAF]=gen_spectrum(x1,x2,step);
I488=imread('新的仿真数据\目标1.jpg');
I555=imread('新的仿真数据\目标2.jpg');
IAF=imread('新的仿真数据\背景.jpg');
chls = size(s488,2);

AFF = 5/10;
d488=double(I488)*10;
d555=double(I555)*10;
dAF=double(IAF)*10*AFF;
% max(max([d488 d555 dAF]))
% figure; imagesc(d488);axis image;axis off;
% figure; imagesc(d555);axis image;axis off;
% figure; imagesc(dAF);axis image;axis off;
% imagesc([d488 d555 dAF]);

[Row Col] = size(d488);
XS = zeros(Row*Col, chls); %用于分离的原始数据
for i = 1:chls
    temp = d488*s488(i)+d555*s555(i)+dAF*sAF(i);
    XS(:,i)=reshape(temp, Row*Col,1);
%     添加噪声
%     temp = d488*s488(i)+d555*s555(i)+dAF*sAF(i);
%     temp = reshape(temp,Row*Col,1);
%     temp = awgn(temp,snr,'measured');
    if i == 3
        h = figure;imagesc(reshape(temp,Row,Col));
        dirname = ['O:\开题报告和学位论文\论文数据\NCut初始化\稀疏系数的影响稀疏通道\orig555_' num2str(spa) 'spa'];
        saveas(h,dirname,'png');
    end
%     XS(:,i)=reshape(temp, Row*Col,1);
end


d488 = imresize(d488,0.25);
d555 = imresize(d555,0.25);
dAF = imresize(dAF,0.25);
[row col] = size(d488);


mask = false(row,col);
mask(dAF > 10) = true;
mask = imopen(mask,strel('disk',1));
% figure,imshow(mask);
m = sum(mask(:));

X = zeros(m,chls);
D = zeros(row,col,chls);
for i = 1:chls
%     temp = d488*s488(i)+d555*s555(i)+dAF*sAF(i);
%     D(:,:,i) = temp;
%     X(:,i) = temp(mask);
   temp = d488*s488(i)+d555*s555(i)+dAF*sAF(i);
%     temp = reshape(temp,row*col,1);
%     temp = awgn(temp,snr,'measured');
%     temp = reshape(temp,row,col);
    D(:,:,i) = temp;
    X(:,i) = temp(mask);
end


%% Normalized Cut进行分割
orignalmask = mask;                    % 保存原始的整个老鼠的掩膜
[y,v,nv,evec,eval] = Ncut(X);
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
eval(2,2)
while eval(2,2) < NcutValue

%     if sum(binpar(:)) > sum(mask(:))/2
%         mask = mask & binpar;
%     else
%         mask = mask & ~binpar;
%     end
    m = sum(mask(:)); 
    X = zeros(m,chls);
    for i = 1:chls
        temp = D(:,:,i);
        X(:,i) = temp(mask);
    end
    [y,v,nv,evec,eval] = Ncut(X);
        eval(2,2)
    if eval(2,2) < NcutValue
        res = zeros(row,col);
        res(mask) = v;
%         figure,imagesc(res);
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

h = figure;imagesc(seg);
dirname = ['O:\开题报告和学位论文\论文数据\NCut初始化\稀疏系数的影响稀疏通道\' num2str(spa) 'spa'];
saveas(h,dirname,'png');
axis image;axis off;
% figure,imshow(mask);

%% 端元提取,其这mask此时此刻已经是纯自发荧光的mask了,targetmask为目标的mask.然后将目标区域减去背景的平均值然后再运用ATGP进行端元提取
targetmask = xor(orignalmask , mask);
figure,imshow(targetmask);
TX = zeros(sum(targetmask(:)), chls);
for i = 1:chls
    temp = D(:,:,i);
    TX(:,i) = temp(targetmask);
end

AverageAuto = zeros(1,chls);
for i = 1:chls
    temp = D(:,:,i);
    AverageAuto(i) = mean(temp(mask));
end
% AverageAuto = [mean(d1(mask)) mean(d2(mask)) mean(d3(mask)) mean(d4(mask)) mean(d5(mask))];
TX = TX - repmat(AverageAuto,sum(targetmask(:)),1);
TX(TX <= 0) = 1e-9;

% [E,C] = EIA_ATGP(TX',2);
% E = E./ repmat(max(E),chls,1);
% figure,plot(x1:step:x2,E(:,1),x1:step:x2,E(:,2),x1:step:x2,AverageAuto/max(AverageAuto));

[E,C] = ATGP(TX,2);
E = E./ repmat(max(E),chls,1);
h = figure;plot(x1:step:x2,E(:,1)*100,x1:step:x2,E(:,2)*100,x1:step:x2,AverageAuto/max(AverageAuto)*100,'LineWidth',3);
axis([x1-10 x2+10 0 105]);
dirname = ['O:\开题报告和学位论文\论文数据\NCut初始化\稀疏系数的影响稀疏通道\s' num2str(spa) 'spa'];
saveas(h,dirname,'png');


% [E, indice, Rp] = VCA(TX','Endmembers',2);
% E = E - repmat(min(E),chls,1);
% E = E./ repmat(max(E),chls,1);
% figure,plot(x1:step:x2,E(:,1),x1:step:x2,E(:,2),x1:step:x2,AverageAuto/max(AverageAuto));
% % figure,plot(1:5,st1,1:5,st2,1:5,sa);
% 
% [E, indice, Rp] = VCA(X','Endmembers',3);
% E = E./ repmat(max(E),chls,1);
% figure,plot(x1:step:x2,E(:,1),x1:step:x2,E(:,2),x1:step:x2,E(:,3));
% 
% [E, C] = EIA_ATGP(X',3);
% E = E./ repmat(max(E),chls,1);
% figure,plot(x1:step:x2,E(:,1),x1:step:x2,E(:,2),x1:step:x2,E(:,3));

%% 接下来就是对原始数据进行分解,这个时候的方法很多
% phi = 0.85;
phi = spa/100;
CompNum = 3;
% C0 = rand(Row*Col, CompNum);
S0 = rand(CompNum,chls);
S0(1:CompNum-1,:) = E';
S0(CompNum, :) = AverageAuto/max(AverageAuto);
[C0,S0]=Initial_Pure(XS,S0,10);
 [C,S,time,cost_record,lambda_sq]= mynmf_ghals_L12_Ultimate(XS,C0,S0,phi);
figure,imagesc(reshape(C(:,1),Row,Col));
figure,imagesc(reshape(C(:,2),Row,Col));
figure,imagesc(reshape(C(:,3),Row,Col));
figure,plot(x1:step:x2,S(1,:),x1:step:x2,S(2,:),x1:step:x2,S(3,:));

save(['O:\开题报告和学位论文\论文数据\NCut初始化\稀疏系数的影响稀疏通道\' num2str(spa) '.mat']);
end 
