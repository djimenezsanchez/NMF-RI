function [IM_reshaped, IM_Mask] = toMask(IM)
 
IM_reshaped = double(reshape(IM,[size(IM,1)*size(IM,2),size(IM,3)]));
% Calculates MASK.
IM_Mask = logical(ones(size(IM_reshaped,1),1));

end