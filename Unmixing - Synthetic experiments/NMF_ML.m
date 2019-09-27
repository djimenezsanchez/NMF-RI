% NMF_ML Blind unmixing based on multi-layering and the L1 norm
%
% NMF_ML(Y,m,A0,H0,split_every,max_it,tolA,display,flatten)
%
%  Y  LxM input mixed matrix
%  m  mask used to determine where the image was sampled
%  A0  LxN initial mixing matrix (optional, defaults to exponential matrix)
%  H0  NxM initial unmixed matrix estimation (optional, defaults to Y)
%  split_every  creates a new layer every split_every iterations (optional)
%  max_it  maximum number of iterations (optional)
%  tolA  tolerance of the A matrix (optional)
%  display  1: display during the iterations, 0: no display
%  flatten  1: for 3D images, calculates MIP before unmixing
%           0: use full volume
%
% Requires DipImage (www.dipimage.org)
%
% Based on the algorithm described in :
% 
% Pengo et al., EFFICIENT BLIND SPECTRAL UNMIXING OF FLUORESCENTLY LABELED 
% SAMPLES USING MULTI-LAYER NON-NEGATIVE MATRIX FACTORIZATION, PLoS ONE, 2013
%
% June 18th 2013  v1.0  First release
% Feb  12th 2014        bugfix affecting flattening, non-rectangular, 3D 
%
function [Ac H xts saveVAR,n] = NMF_ML(Y,m,varargin)
n=size(Y,1);

p = inputParser;
p.addOptional('A0',gauss_matrix(n), @(A) size(A,1)==size(Y,1));
p.addOptional('H0',Y, @(H) size(H,2)==size(Y,2));
p.addOptional('split_every',20);
p.addOptional('max_it',1500);
p.addOptional('tolA',1e-2);
p.addOptional('display',0);
p.addOptional('flatten',1);
p.addOptional('alphaA',0.1);
p.addOptional('alphaX',0.01);
p.parse(varargin{:});
s = p.Results;

saveVAR=[];

% Initialize
A0=s.A0;
H0=s.H0;
A0 = A0*diag(1./(sum(A0,1)+eps)); % DJ
A = A0;
H=H0;

% Transform arguments to variables
split_every = s.split_every;
max_it = s.max_it;
tolA = s.tolA;
display = s.display;
flatten = s.flatten;

alphaA = s.alphaA;
alphaX = s.alphaX;
psi = 1;

if flatten
    if size(m,3)>1
        % Flatten (calculate MIP)
        im=iterate('max',dip_image(reconstruct_image(m,Y)),[],3);
        im1=iterate('max',dip_image(reconstruct_image(m,H)),[],3);
        m=max(m,[],3);
        Y=[]; H=[];
        for i=1:length(im)
            Y=cat(1,Y,double(im{i}(m)));
        end
	for i=1:length(im1)
            H=cat(1,H,double(im1{i}(m)));
        end
    end
end

Y0=Y;
    
cXT = @(a,B) mean(B(a))/mean(B);
xts = [];
for n=1:max_it
       
    if mod(n,split_every)==0
        if n==split_every
            Ac=A;
            A=exponential_matrix(size(Ac,2));
        else
            Ac=Ac*A;
        end
        
        Y=H;

        H=max(1E6*eps,pinv(Ac'*Ac)*Ac'*Y0);
        
        % Max distance from 0 or 1 < tolA       
        fprintf('M: %02.4f   Fro: %02.4f\n', max(.5-abs(A(:)-.5)), norm(A-eye(size(A)),'fro'))
        if max(.5-abs(A(:)-.5))<tolA
            xt = cXT(Ac(end-1,end-1)*Y(end-1,:)>Ac(end,end-1)*Y(end,:),H(end,:));
            xts=cat(1,xts,xt);
            break
        end
    end
    
    % H step
    Yx = A'*Y - alphaX*psi;
    Yx(Yx <= 0) = 100*eps;
    H = H.*(Yx./((A'*A)*H + eps));

    % A step
    %Ap = A;
    Ya = H*Y' - alphaA;
    Ya(Ya <= 0) = 100*eps;
    A = A.*(Ya./((A*(H*H'))' + eps))';
    A = A*diag(1./(sum(A,1)+eps));

    if ~exist('Ac','var')
        Ac=A;
    end

    xt = cXT(Ac(end-1,end-1)*Y(end-1,:)>Ac(end,end-1)*Y(end,:),H(end,:));
    xts=cat(1,xts,xt);
    
    if display
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %Z = A*H + eps;
        %Z = diag(1./sqrt(var(Z,2)'))*Z;
        %kl = sum(sum(Y.*log(Y./Z + eps) - Y + Z));
        if mod(n,20)==0
            subplot(1,4,1)
                plot(Y(end-1,:),Y(end,:),'.',H(end-1,:),H(end,:),'k.')
            subplot(1,4,2)
                plot(xts)
                text(n/4,mean(xts),sprintf('%5d %1.5f',n,xt))
            if ~isempty(m)
                subplot(1,4,3)
                im=reconstruct_image(m,H);
                M=zeros(size(im{1},2),size(im{1},1),3);
                for i=1:length(im)
                    if size(m,3)>1
                        mp = max(im{i},[],3);
                    else
                        mp = m;
                    end
                    if i==4
                        for j=1:3
                            M(:,:,j)=min(1,M(:,:,j)+double(mp/max(mp)));
                        end
                    else
                        M(:,:,i)=double(mp/max(mp));
                    end
                end
                image(M)
            end
            subplot(1,4,4)
                if exist('Ac','var')
                    plot(Ac)
                else
                    plot(A)
                end
            drawnow
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end

H=max(1E6*eps,pinv(Ac'*Ac)*Ac'*Y0);

function im=reconstruct_image(m,Y)

Nch = size(Y,1);
im = newimar(Nch);
for i=1:Nch
    im{i}=newim(size(m));
    im{i}(m)=Y(i,:);
end


