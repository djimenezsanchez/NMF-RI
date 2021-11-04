% NMF_RI Blind unmixing based on Reinitialization
%
% NMF_RI(Y,A0,H0)
%
%  Y  LxM input mixed matrix
%  A  LxN initial mixing matrix
%  H  NxM initial unmixed matrix estimation
%  max_it  maximum number of iterations
%
% Parts of this code were forked from T.Pengo code described in: 
% "Efficient Blind Spectral Unmixing of Fluorescently Labeled Samples Using
% Multi-Layer Non-Negative Matrix Factorization"
% 
% This algorithm is described in:  
% "NMF-RI: Blind spectral unmixing of highly mixed multispectral flow and image cytometry"


function [A H] = NMF_RI(Y,A0,H0)

% Initialize
A0 = A0*diag(1./(sum(A0,1)+eps));
H=H0;

% Variables
max_it = 300; % Number of maximum iterations
theta1 = 0.01; % Stopping Criteria intraLayer
theta2 = 0.01; % Stopping Criteria interLayer
alphaA = 0.001; % Sparseness for A
alphaX = mean(Y(:))*0.001; % Sparseness for H
    
layer=1;
LayerN = 0;

Astep=A0; % between steps
Aprev=A0; % inside step
A = A0;
for n=1:max_it
    
    Anow = A; 
    
    aVect(n-LayerN) = sum(sum(abs(Aprev-Anow)));
    if n==1
        aVect(1)=100;
    end
    daVect = diff(aVect);
     
    if n>1
    if any(abs(aVect)<theta1) | any(daVect>0) %Optimization triggering                     
        disp(['converge: ',num2str(sum(sum(abs(Astep-A)))),' N: ',num2str(n)]); 

        if layer==1            
            Astep=A;
            A = A*diag(1./(sum(A,1)+eps));   
        else             
            if sum(sum(abs(Astep-A)))<theta2
            break
            else
                Astep=A;
                A = A*diag(1./(sum(A,1)+eps));
            end            
        end
        
        layer=layer+1;  
        LayerN = n;
        clear aVect; 
         
        H=max(1E6*eps,pinv(Astep'*Astep)*Astep'*Y);        
    end
    end    
    Aprev = A;
        
    % H step    
    Yx = A'*Y - alphaX;
    H = H.*(Yx./((A'*A)*H + eps));
    
    % A step
    Ya = H*Y' - alphaA;
    A = A.*(Ya./((A*(H*H'))' + eps))';
%     A(:,end) = A0(:,end); %Fix AF

    A = A*diag(1./(sum(A,1)+eps));    
end


H=max(1E6*eps,pinv(A'*A)*A'*Y);

end
