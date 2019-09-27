% Data preprocessing method as described in Algorithm 1. 

function [Y_sps, H_sps] = dataPreprocessing(A0,Y)

% Calculate H0 using A and Y.
H0 = max(1E6*eps,pinv(A0'*A0)*A0'*Y);

% Choose Sparseness selection method or most similar events to each
% spectra.
sparsenessMethod=false; 
Ths = 0.99; % Take 1% of sparsest pixels
if sparsenessMethod
    SPS = ( sqrt(size(H0,1))-(sum(abs(H0),1)./sqrt(sum(H0.^2,1))) )/(sqrt(size(H0,1))-1);                
    [f, x] = ecdf(SPS);
    thresh = x(round(length(x)*Ths));
    Y_sps=Y(:,SPS>thresh);
    Y_sps = Y_sps./mean(Y_sps,1);
    H_sps = max(1E6*eps,pinv(A0'*A0)*A0'*Y_sps);

else
    % Similarity measure between signal and each theoreticalSpectra
    logicalMap = logical(zeros(1,size(Y,2)));
    % normalize using mean to allow comparison between cells.
    Y_to_Segment = Y./mean(Y,1);
    AtheoricalSim = A0./mean(A0); 

    % Take the most similar events to each spectra. 
    for i = 1:size(A0,2)-1 % All fluorochromes
        % similarity measure 
        RMSE(:,i) = sum((Y_to_Segment.*repmat(AtheoricalSim(:,i),[1,size(Y,2)])),1);  % Root Mean Squared Error bet
        % Threshold values
        [counts,x] = imhist(RMSE(:,i)./max(RMSE(:,i)),100);
        % save logical values
        [f,x] = ecdf(RMSE(:,i)./max(RMSE(:,i)));
        thresh = x(round(length(x)*Ths)); 
        logicalMap = logicalMap | logical(thresh<RMSE(:,i)./max(RMSE(:,i)))';
    end
    Y = Y(:,logicalMap); % Selected pixels.
    Y_sps = Y./mean(Y,1);
    H_sps = max(1E6*eps,pinv(A0'*A0)*A0'*Y_sps);
end
end