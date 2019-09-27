clear d488 d555 dAF

I488=imread('新的仿真数据\目标1.jpg');
I555=imread('新的仿真数据\目标2.jpg');
IAF=imread('新的仿真数据\背景.jpg');
d488=double(I488)*10*max(s488)/100;
d555=double(I555)*10*max(s555)/100;
dAF=double(IAF)*10*AFF*max(sAF)/100;


% RMSE_3D([reshape(d488/max(d488(:)),Row*Col,1)...
%     reshape(d555/max(d555(:)),Row*Col,1)...
%     reshape(dAF/max(dAF(:)), Row*Col,1)], [C(:,2)/max(C(:,2)) C(:,1)/max(C(:,1)) C(:,3)/max(C(:,3))])
%   
An = [reshape(d555,Row*Col,1) reshape(d488,Row*Col,1) reshape(dAF,Row*Col,1)];
CC = C.*(ones(Row*Col,1)*max(S,[],2)');
CC=CC/100;
RMSE_3D(An,CC)