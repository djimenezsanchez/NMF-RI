%  保存图像
for j = 60:5:85
    load(['O:\开题报告和学位论文\论文数据\NCut初始化\稀疏系数的影响\' num2str(j) '.mat']);
    dirname=['O:\开题报告和学位论文\论文数据\NCut初始化\稀疏系数的影响\' num2str(j)];%新的文件夹名
    a=['mkdir ' dirname];%创建命令
    system(a) %创建文件夹
    CC = C.*(ones(Row*Col,1)*max(S,[],2)');
%     CC=CC/100;
%     V = max(C(:));
    V = 2800;
    h = figure;imagesc(reshape(CC(:,1),Row,Col), [0 V]);axis off;axis image;
%     h = figure;image(reshape(CC(:,1)/V,Row,Col));axis off;axis image;
%     colormap(jet);
    saveas(h,[dirname '\af488'],'png');
    disp('第1个稀疏系数');
    Sparsity_Vec(CC(:,1),0)
    disp('\n');
    h = figure;imagesc(reshape(CC(:,2),Row,Col),[0 V]);axis off;axis image;
%     h = figure;image(reshape(CC(:,2),Row,Col)/V);axis off;axis image;
%     colormap(jet);
    saveas(h,[dirname '\af555'],'png');
    disp('第2个稀疏系数');
    Sparsity_Vec(CC(:,2),0)
    disp('\n');
    h = figure;imagesc(reshape(CC(:,3),Row,Col),[0 V]);axis off;axis image;
%     h = figure;image(reshape(CC(:,3)/V,Row,Col));axis off;axis image;
%     colormap(jet);
    saveas(h,[dirname '\af'],'png');
    disp('第3个稀疏系数');
    Sparsity_Vec(CC(:,3),0)
    disp('\n');
    W = repmat(max(S,[],2),1,chls);
    S = S./W;
    h = figure;plot(x1:step:x2,S(1,:)*100,x1:step:x2,S(2,:)*100,x1:step:x2,S(3,:)*100,'LineWidth',3);
    axis([470 660 0 105]);
    xlabel('Wavelength(nm)')
    legend('AF555','AF488','Autofluorescence');
    saveas(h,[dirname '\spectrum'],'png');
    close all;
    clear all;
end 

% load(['O:\开题报告和学位论文\论文数据\NCut初始化\HALSL1的结果数据\halsL0_01.mat']);
% dirname=['O:\开题报告和学位论文\论文数据\NCut初始化\HALSL1的结果数据\'];%新的文件夹名
% a=['mkdir ' dirname];%创建命令
% system(a) %创建文件夹
% CC = C.*(ones(Row*Col,1)*max(S,[],2)');
% %     CC=CC/100;
% %     V = max(C(:));
% V = 3200;
% h = figure;imagesc(reshape(CC(:,1),Row,Col), [0 V]);axis off;axis image;
% %     h = figure;image(reshape(CC(:,1)/V,Row,Col));axis off;axis image;
% %     colormap(jet);
% saveas(h,[dirname '
% \af488'],'png');
% disp('第1个稀疏系数');
% Sparsity_Vec(CC(:,1),0)
% disp('\n');
% h = figure;imagesc(reshape(CC(:,2),Row,Col),[0 V]);axis off;axis image;
% %     h = figure;image(reshape(CC(:,2),Row,Col)/V);axis off;axis image;
% %     colormap(jet);
% saveas(h,[dirname '\af555'],'png');
% disp('第2个稀疏系数');
% Sparsity_Vec(CC(:,2),0)
% disp('\n');
% h = figure;imagesc(reshape(CC(:,3),Row,Col),[0 V]);axis off;axis image;
% %     h = figure;image(reshape(CC(:,3)/V,Row,Col));axis off;axis image;
% %     colormap(jet);
% saveas(h,[dirname '\af'],'png');
% disp('第3个稀疏系数');
% Sparsity_Vec(CC(:,3),0)
% disp('\n');
% W = repmat(max(S,[],2),1,chls);
% S = S./W;
% h = figure;plot(x1:step:x2,S(1,:)*100,x1:step:x2,S(2,:)*100,x1:step:x2,S(3,:)*100,'LineWidth',3);
% axis([470 660 0 105]);
% xlabel('Wavelength(nm)')
% legend('AF555','AF488','Autofluorescence');
% saveas(h,[dirname '\spectrum'],'png');
% % close all;
% % clear all;
% 
% 
