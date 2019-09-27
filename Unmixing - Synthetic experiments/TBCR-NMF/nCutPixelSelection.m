function [] = nCutPixelSelection(Y)

tic;
orignalmask = mask;                    % the mask of whole mouse
[y,v,nv,evec,eval] = Ncut(Y);%nv
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
while eval(2,2) < NcutValue
%     if sum(binpar(:)) > sum(mask(:))/2
%         mask = mask & binpar;
%     else
%         mask = mask & ~binpar;
%     end
    eval(2,2)
    m = sum(mask(:)); 
    Y = zeros(m,4);
    Y(:,1) = d1(mask);
    Y(:,2) = d2(mask);
    Y(:,3) = d3(mask);    

    [y,v,nv,evec,eval] = Ncut(Y);%nv
    if eval(2,2) < NcutValue
        res = zeros(row,col);
        res(mask) = v;
        figure,imagesc(res);
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

