function [y,v,nv,evec,eval] = Ncut(X)

[m,n] = size(X);
%每一行中心化，归一化
% CX = X - repmat(mean(X,2),1,n); 
% CX = CX./ repmat(sqrt(sum(CX.^2,2)),1,n);
CX = X./(max(X,[],2)*ones(1,n));

sums = sum(CX.*CX,2);
K = sums * ones(1,m);
W = K + K' - 2*(CX*CX');
% sigma_I = 0.06;% 仿真数据这个值是这个
% sigma_I = 1;
sigma_I = 0.1;
% sigma_I = 0.08;
W = exp(-W/sigma_I^2);
% W = zeros(m,m);
% sigma_I = 0.01;
% for i = 1:m
%     for j = 1:m
% %         distance = sqrt((CoordX(i) - CoordX(j))^2 + (CoordY(i) - CoordY(j))^2);
% %         if distance < r
% %             W(i,j) = exp(- (CX(i,:) - CX(j,:))*(CX(i,:) - CX(j,:))'/sigma_I);
% %             W(i,j) = W(i,j) * exp(-distance^2/sigma_X);
% %         else
% %             W(i,j) = 0;
% %         end
%     W(i,j) = exp(- (CX(i,:) - CX(j,:))*(CX(i,:) - CX(j,:))'/sigma_I);
% %         W(i,j) = exp(- (DX(i) - DX(j))*(DX(i) - DX(j))'/sigma_I);
%     end
% end
W = (W + W')/2;
D = diag(sum(W));
% figure,imagesc(W);
% 计算特征向量
% [evec eval] = eig(D-W,D);

T = D-W;
T = (T+T')/2;
opt = struct('issym',true,'isreal',true);
[evec eval] = eigs(T,D, 2, 'sa', opt );
% [evec eval] = eig(T,D);
v = evec(:,2);

k = sum(sum(W(v>0,:)))/sum(sum(W));
y = zeros(m,1);
y(v>0) = 1;
y(v<=0) = -k/(1-k);
%  y(v<=0)=-1;
nv = y'*(D-W)*y/(y'*D*y);
end


