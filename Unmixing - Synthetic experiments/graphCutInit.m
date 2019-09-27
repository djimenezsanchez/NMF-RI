function [E] = graphCutInit(Y,numFluorochromes)
% Forked from TBCR-NMF
% Qin, Binjie & Hu, Chen & Huang, Shaosen. (2016).
% Target/Background Classification Regularized Nonnegative Matrix Factorization for Fluorescence Unmixing.
% IEEE Transactions on Instrumentation and Measurement. 65. 1-16. 10.1109/TIM.2016.2516318. 
TBCRY_raw = Y(:,1:1:end)';
mask = logical(ones(size(TBCRY_raw,1),1));
reg = 0;
seg = zeros(size(TBCRY_raw,1),1);
for i = 1:numFluorochromes
    TBCRY = TBCRY_raw(mask,:);
    [y,v,nv,evec,eval] = Ncut(TBCRY);%nv
    res = zeros(size(mask));
    res(mask) = v;
    binpar = true(size(mask));
    binpar(res<0) = false;
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
[E,C] = EIA_ATGP(TBCRY_raw(seg>0,:)',numFluorochromes);    
[~,maxI]=max(E);
[~,howtoSort] = sort(maxI);
E = E(:,howtoSort);


