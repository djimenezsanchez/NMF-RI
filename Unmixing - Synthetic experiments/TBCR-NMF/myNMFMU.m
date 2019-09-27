function [C,S]=myNMFMU(D,C,S,maxiter)
% S=Pure;
% C=rand(size(D,1),size(Pure,1));
% for i=1:100 %% 2_26之前的版本是20
%     C=C.*(D*S')./(C*(S*S' )+ 1e-9);
% end;
for i=1:maxiter %% 2_26之前的版本是20
    C=C.*(D*S')./(C*(S*S' )+ 1e-9);
    S=S.*(C'*D)./(C'*C*S+1e-9);
end;
