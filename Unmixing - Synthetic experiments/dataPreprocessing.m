function [Y,Y_sps, H_sps] = dataPreprocessing(A0,Y)

% Elimination of dark current.
Y = Y-min(Y);

% Matrix Initialization
H0 = max(100*eps,pinv([A0])*[Y]); 

% Sparseness selection method.
Thrs=0.7; % Sparseness threshold, select 30% of pixels.
H_sps = H0;
Y_sps = Y;
% segmentation to obtain sparsest pixels
SPS = ( sqrt(size(H_sps,1))-(sum(abs(H_sps),1)./sqrt(sum(H_sps.^2,1))) )/(sqrt(size(H_sps,1))-1);                
% Cumulative distribution function to find sparsest pixels.
[f, x] = ecdf(SPS);
thresh = x(round(length(x)*Thrs));
Y_sps=Y(:,SPS>thresh);  
% Prepare H0 to unmix.-
H_sps = max(100*eps,pinv(A0)*Y_sps); 

end