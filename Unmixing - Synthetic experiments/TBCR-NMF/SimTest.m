close all;clear all;
x1 = 480;x2=650;step=5;
[s488 s555 sAF]=gen_spectrum(x1,x2,step);
I488=imread('simulation\target1.jpg');
I555=imread('simulation\target2.jpg');
IAF=imread('simulation\background.jpg');
chls = size(s488,2);

AFF = 3/10;
d488=double(I488)*10;
d555=double(I555)*10;
dAF=double(IAF)*10*AFF;
% max(max([d488 d555 dAF]))
% figure; imagesc(d488);axis image;axis off;
% figure; imagesc(d555);axis image;axis off;
% figure; imagesc(dAF);axis image;axis off;
% imagesc([d488 d555 dAF]);

[Row Col] = size(d488);
XS = zeros(Row*Col, chls); % the original source data
for i = 1:chls
    temp = d488*s488(i)+d555*s555(i)+dAF*sAF(i);
    XS(:,i)=reshape(temp, Row*Col,1);
%     add noise
% for snr = 15:5:35
%     temp = d488*s488(i)+d555*s555(i)+dAF*sAF(i);
%     temp = reshape(temp,Row*Col,1);
%     temp = awgn(temp,snr,'measured');
    if i == 15
        h = figure;imagesc(reshape(temp,Row,Col),[1 2550]);axis off;axis image;
%         dirname = ['O:\SNREffect\orig555_' num2str(snr) 'db'];
%         saveas(h,dirname,'png');
    end
    XS(:,i)=reshape(temp, Row*Col,1);
% end
end
figure,imagesc(d488);axis image;axis off;
figure,imagesc(d555);axis image;axis off;
figure,imagesc(dAF);axis image;axis off;

d488 = imresize(d488,0.25);
d555 = imresize(d555,0.25);
dAF = imresize(dAF,0.25);
[row col] = size(d488);


mask = false(row,col);
mask(dAF > 10) = true;
% mask = imopen(mask,strel('disk',1));
mask = imerode(mask,strel('disk',1));
% figure,imshow(mask);
m = sum(mask(:));

X = zeros(m,chls);
D = zeros(row,col,chls);
for i = 1:chls
    temp = d488*s488(i)+d555*s555(i)+dAF*sAF(i);
    D(:,:,i) = temp;
    X(:,i) = temp(mask);
%    temp = d488*s488(i)+d555*s555(i)+dAF*sAF(i);
%     temp = reshape(temp,row*col,1);
%     temp = awgn(temp,snr,'measured');
%     temp = reshape(temp,row,col);
%     D(:,:,i) = temp;
%     X(:,i) = temp(mask);
end


%% Normalized Cut for target-background classification
orignalmask = mask;                    % save the mask region of whole mouse
[y,v,nv,evec,eval] = Ncut(X);
res = zeros(row,col);
res(mask) = v;
% figure,imagesc(res);

binpar = false(row,col);
binpar(res>0.0) = true;
% figure,imshow(binpar);
binpar = imopen(binpar,strel('square',1));

seg = zeros(row,col);
reg = 1;
seg(mask) = reg;
reg = reg + 1;
if sum(binpar(:)) > sum(mask(:))/2
    seg(mask & ~binpar) = reg;
    mask = mask & binpar;
    mask = imopen(mask,strel('square',1));
else
    seg(mask & binpar) = reg;
    mask = mask & ~binpar;
    mask = imopen(mask,strel('square',1));
end

NcutValue = 0.4;
eval(2,2)
count = 1;
while eval(2,2) < NcutValue && count < 4
    count = count +1;
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
% dirname = ['O:\SNREffect\' num2str(snr) 'db'];
% saveas(h,dirname,'png');
axis image;axis off;
% figure,imshow(mask);

%% Endmember extraction using ATGP after the average background autofluorescence are removed from the target regions
%% the mask is the mask from the background autofluorescence¡À,targetmask is the mask for the target regions.
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
h = figure;
plot(x1:step:x2,E(:,1)*100,'color',[0,0.5,0],'LineWidth',2);
hold on;plot(x1:step:x2,E(:,2)*100,'b','LineWidth',2);
hold on;plot(x1:step:x2,AverageAuto/max(AverageAuto)*100,'r','LineWidth',2);
legend('AF488','AF555','Autofluorescence');

axis([x1-10 x2+10 0 105]);

% dirname = ['O:\SNREffect\s' num2str(snr) 'db'];
% saveas(h,dirname,'png');


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

%% unmixing the original mixed data
phi = 0.75; CompNum = 3;
% C0 = rand(Row*Col, CompNum);
S0 = rand(CompNum,chls);
% S0(1:CompNum-1,:) = E';
% S0(CompNum, :) = AverageAuto/max(AverageAuto);
[C0,S0]=Initial_Pure(XS,S0,100);
%  [C,S,time,cost_record,lambda_sq]= mynmf_ghals_L12_Ultimate(XS,C0,S0,phi);
% [C,S,time,cost_record]= mynmf_ghals_L1(XS,C0,S0,phi);
% [C,S,t,objhistory]=myNMFsv(XS,phi-0.1,C0,S0);%,tol)


% save(['O:\SNREffect\' num2str(snr) '.mat']);
% [C,S,conv] = hc_sparseNMF(XS,phi,C0,S0,2);
% C0 = rand(Row*Col,3);
% S0 = rand(3,chls);
[C,S,time,cost_record]= mynmf_ghals(XS,C0,S0,20);
% [C,S,t,objhistory]=myNMFsv(XS,0.85,C0,S0);
%  [C,S,time,cost_record]= mynmf_ghals_L1(XS,C0,S0,0.85);
% [C,S] = myHALS(XS,C0,S0);
%  [C,S,time,cost_record,lambda_sq] = mynmf_ghals_L12_Ultimate(XS,C0,S0,phi);
% S=S./repmat(sqrt(sum(S.^2,2)),[1 size(S,2)]);
% C=C.*repmat(sqrt(sum(S.^2,2))',[size(C,1) 1]);
C = C.*(repmat(max(S,[],2)',Row*Col,1));
S = S./(max(S,[],2)*ones(1,chls));
figure,imagesc(reshape(C(:,1),Row,Col));axis off;axis image;
figure,imagesc(reshape(C(:,2),Row,Col));axis off;axis image;
figure,imagesc(reshape(C(:,3),Row,Col));axis off;axis image;
% figure,plot(x1:step:x2,S(1,:),x1:step:x2,S(2,:),x1:step:x2,S(3,:));
% axis([480 675 0 1.1]);
lw = 2;
figure,plot(x1:step:x2, S(1,:)*100,'b','LineWidth',lw);
hold on;plot(x1:step:x2, S(2,:)*100,'c','LineWidth',lw, 'Color',[0,0.5,0]);
hold on;plot(x1:step:x2, S(3,:)*100,'r','LineWidth',lw);
axis([470 660 0 105]);
legend('AF555','AF488','Autofluorescence');
xlabel('Wavelength(nm)');
ylabel('Relative Intensity');
% close all