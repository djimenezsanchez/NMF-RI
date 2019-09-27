% test spectral cluster
x = 0:0.1:6.2;
x = x + 1*rand(1,size(x,2));

r1 = 4;

x1 = (r1)*cos(x);
y1 =  (r1)*sin(x);


r2 = 8;
x2 = r2*cos(x);
y2 = r2*sin(x);

X = [x1 x2;y1 y2];
plot(X(1,:),X(2,:),'g*');

X = X';

m = size(X,1);

W = zeros(m,m);
r = 5;sigma_I = 2;
for i = 1:m
    for j = 1:m
%         distance = sqrt((CoordX(i) - CoordX(j))^2 + (CoordY(i) - CoordY(j))^2);
%         if distance < r
%             W(i,j) = exp(- (CX(i,:) - CX(j,:))*(CX(i,:) - CX(j,:))'/sigma_I);
%             W(i,j) = W(i,j) * exp(-distance^2/sigma_X);
%         else
%             W(i,j) = 0;
%         end
    W(i,j) = exp(- (X(i,:) - X(j,:))*(X(i,:) - X(j,:))'/sigma_I);
%         W(i,j) = exp(- (DX(i) - DX(j))*(DX(i) - DX(j))'/sigma_I);
    end
end
D = diag(sum(W));
figure,imshow(W,[]);
% 计算特征向量
[evec eval] = eig(D-W,D);
v = evec(:,2);

idx = find(v > 0);
figure,
plot(X(idx,1), X(idx,2),'go');
hold on
idx = find(v < 0);
plot(X(idx,1), X(idx,2),'ro');