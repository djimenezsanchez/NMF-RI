function [U,index]=ATGP(data,p)
tic;
%X=data';
X=data;
[r,c]=size(X);%rÎªN£¬cÎªL
I=eye(c);
XX=X.*X;
XX_row=sum(XX,2);
[XX_max,index]=sort(XX_row,'descend');
k=index(1,1);
t0=X(k,:)';
X(k,:)=0;
clear XX XX_row index XX_max
P_t=I-t0*inv(t0'*t0)*t0';
for i=1:p
    PX=X*P_t';
    PXPX=PX.*PX;
    PXPX_row=sum(PXPX,2);
    [PXPX_max,index1]=sort(PXPX_row,'descend');
    k=index1(1,1);
    index(i,1)=k;
    U(:,i)=X(k,:)';
    X(k,:)=0;
    P_t=I-U*inv(U'*U)*U';
end;
clear P_t PX PXPX 
toc
end