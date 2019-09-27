function [C0,S0]=Initial_Pure(D,Pure,maxiter)
S0=Pure;
C0=rand(size(D,1),size(Pure,1));
for i=1:maxiter %% 2_26之前的版本是20
    C0=C0.*(D*S0')./(C0*(S0*S0' )+ 1e-9);
end;