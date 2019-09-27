%  Lagrangian NMU algorithm
%
% cf. "Using Underapproximations for Sparse Nonnegative Matrix Factorization",
% N. Gillis and F. Glineur, Pattern Recognition 43 (4), p.1676-1687, 2010.
% website. http://sites.google.com/site/nicolasgillis/
%
%
% [V,W,e] = LagrangianNMU(M,r,maxiter,Vinit,Winit)
%
% Input.
%   M              : (m x n) matrix to factorize
%   maxiter        : number of iterations
%   r              : factorization rank
%   (Vinit, Winit) : initial matrices of dimensions (m x r) and (r x n)
%
% Output.
%   (V,W) : solution, VW "<=" M, V >= 0, W >= 0   ("." - approximately)
%   e     : evolution of ||M-VW||_F^2

function [V,W,e] = LagrangianNMU(M,r,maxiter,Vinit,Winit)

[m,n] = size(M);
% Initialization
if nargin == 3
    V = rand(m,r); W = rand(r,n);
else
    V = Vinit; W = Winit;
end
% Scaling
A = V*W;
alpha = sum(sum(A.*M))/norm(A,'fro')^2;
V = V*alpha;
% Initialization of Lagrangian variable to 0
lambda = zeros(m,n);
% Parameter rho, can be tuned to improve performances
rho0 = 1; rho = rho0;

% Main loop
for j = 1 : maxiter
    % ***HALS Update***
    for p = 1 : 2  % Parameter: number of HALS updates between each update of Lambda.
        % Update V
        Ml = (M-lambda);
        A = Ml*W'; B = W*W';
        for k = 1 : r
            V(:,k) = max(A(:,k)-V(:,[1:k-1 k+1:end])*B([1:k-1 k+1:end],k),0)/(B(k,k)+1e-16);
            % If V(:,k) == 0 : numerical problems and rank deficient
            % solution. If it happens, set the vector to a small constant. It
            % works well in practice.
            if V(:,k) == 0
                V(:,k) = 1e-16*ones(m,1);
            end
        end
        % Update W
        A = V'*Ml; B = V'*V;
        for k = 1 : r
            W(k,:) = max(A(k,:)-B(k,[1:k-1 k+1:end]) * W([1:k-1 k+1:end],:),0) /(B(k,k)+1e-16);
            if W(k,:) == 0
                W(k,:) = 1e-16*ones(1,n);
            end
        end
    end
    %   % *** For MU updates ***, use
    %   for p = 1 : 2
    %       W = W.*(V'*(M))./((V'*V*W+V'*lambda)+1e-16);
    %       V = V.*((M)*W')./( (V*(W*W')+lambda*W')+1e-16);
    %   end

    % normalization : ||V_{:k}|| = ||W_{k:}||, forall k.
    norm2v=sqrt(sum(V.^2,1)); norm2w=sqrt(sum((W').^2,1));
    norm2vw = sqrt(norm2v.*norm2w);
    V=V.*repmat(norm2vw./norm2v,m,1);  W=W.*repmat(norm2vw'./norm2w',1,n);
    % ***Lambda Update***
    R = M-V*W; % residue
    lambda = max(0,lambda - rho*R);
    % Error
    e(j) = norm(R,'fro')^2;
    % Update of rho
    rho = rho0/j;
end