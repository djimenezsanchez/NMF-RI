function [W, H, errorHistory] = sparseNMF_W(X, K, L, numIter, numNMFIter, DISPLAY)
% [W, H, errorHistory] = sparseNMF_W(X, K, L, numIter, numNMFIter, DISPLAY)
%
% Sparse NMF with ℓ0-constraints on dictionary matrix W.
% Minimizes approximately sum(sum((X - W*H).^2)) w.r.t. W and H,
% s.t. all(H(:) >= 0), all(W(:) >= 0) and all(sum(W > 0,1) <= L). 
% I.e. X is approximated by W*H, where maximal L entries are nonzero
% in each column of W.
% X has to be nonnegative, i.e. all(X(:) >= 0).
%
% Input:
% X: DxN nonnegative data matrix
% K: inner approximation rank (W -> DxK, H -> KxN)
% L: sparseness factor (maximal number of nonzeros per column of W)
% numIter: number of outer iterations
% numNMFIter: number of NMF iterations for coding matrix update stage
% DISPLAY: when ~= 0, progress is displayed
%
% Output:
% W: DxK sparse nonnegative dictionary
% H: KxN coding matrix
% errorHistory: error as function of iteration in terms of RMSE
%
%
% see "Sparse Nonnegative Matrix Factorization Using ℓ0-constraints",
% R. Peharz, M. Stark, F. Pernkopf, MLSP 2010.
%
% March 2010, by Robert Peharz


[D,N] = size(X);
errorHistory = zeros(numIter,1);

H = rand(K,N);

for k=1:numIter
    % this approximates nonnegative least squares (NMF for W)
    W = ones(D,K);
    for i=1:numNMFIter
        W = W .* (X*H') ./ (W*H*H' + 1e-9);
    end
    
    % Set K-L smallest values to zero for each atom
    for p=1:K
        [sw, idx] = sort(W(:,p),'ascend');
        W(idx(1:D-L),p) = 0;
    end
    
    % NMF
    for i=1:numNMFIter
        W = W .* (X*H') ./ (W*H*H' + 1e-9);
        H = H .* (W'*X) ./ (W'*W*H + 1e-9);
    end
    
    errorHistory(numIter) = sqrt(sum(sum((X-W*H).^2)) / (D*N));
    
    if DISPLAY
        fprintf('Iteration: %d   Error: %f\n',k, errorHistory(end));
    end
end


