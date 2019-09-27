% Synthetic test and performance comparison between algorithms (NMF-RI, NMF-ML and NMF-SB)
% It follows NMF-RI as described in the paper "NMF-RI: Blind spectral unmixing of highly mixed multispectral flow and image cytometry"
% Observation matrices were obtained from Vectra Polaris as a tiff-stack.
% We used emission and excitation data to create theoretical spectra.
% Written by Daniel Jimenez-Sanchez.

clear all;
close all; 
clc;

% To execute different studies:
% Choose between -> 'Colocalization', 'RMSEThRef', 'SNR', 'SpectralResolution'
Study = 'RMSEThRef'; 
    
switch Study
    case 'Colocalization'
        %%%%%%%%%%%%%%%%%%%%%%%-Parameters-%%%%%%%%%%%%%%%%%%%%%%%
        % These parameters studies different levels of colocalization [0,25,50,75 and 100].
        nf = 3;                           % number of fluorochromes        
        VectSNR = 35;                     % SNR levels
        VectRMSEThRef = 3.5;              % Similarity between Reference and theoretical spectra
        Vectnfiles =1:5;                 % Files to open with different levels of colocalization
        VectnChannels=36;                 % Number of channels              
        ChoosingMethod = [1,2,3,4,5,6,13,14,16]; % Number of methods to compare        
        %%%%%%%%%%%%%%%%%%%%%%%-Parameters-%%%%%%%%%%%%%%%%%%%%%%%
    case 'RMSEThRef'
        %%%%%%%%%%%%%%%%%%%%%%%-Parameters-%%%%%%%%%%%%%%%%%%%%%%%
        % These parameters studies different RMSE between Theoretical and Reference Spectra.
        nf = 3;                     % number of fluorochromes        
        VectSNR = 35;               % SNR levels
        VectRMSEThRef = [2:0.5:4.5, 4.25];    % Similarity between Reference and theoretical spectra
        Vectnfiles = [4];             % Files to open with different levels of colocalization
        VectnChannels=36;           % Number of channels              
        ChoosingMethod = [1,2,3,5,13]; % Number of methods to compare        
        %%%%%%%%%%%%%%%%%%%%%%%-Parameters-%%%%%%%%%%%%%%%%%%%%%%%
    case 'SNR'
        %%%%%%%%%%%%%%%%%%%%%%%-Parameters-%%%%%%%%%%%%%%%%%%%%%%%
        % These parameters studies different levels of SNR.
        nf = 3;                     % number of fluorochromes        
        VectSNR = [20:5:40];        % SNR levels
        VectRMSEThRef = 3.5;        % Similarity between Reference and theoretical spectra
        Vectnfiles = [4:5];             % Files to open with different levels of colocalization
        VectnChannels=36;           % Number of channels              
        ChoosingMethod = [1,2,3,5,8,9,10,13,15]; % Number of methods to compare        
        %%%%%%%%%%%%%%%%%%%%%%%-Parameters-%%%%%%%%%%%%%%%%%%%%%%%
    case 'SpectralResolution'
        %%%%%%%%%%%%%%%%%%%%%%%-Parameters-%%%%%%%%%%%%%%%%%%%%%%%
        % These parameters studies different levels of SNR.
        nf = 3;                     % number of fluorochromes        
        VectSNR = 35;               % SNR levels
        VectRMSEThRef = 3.5;        % Similarity between Reference and theoretical spectra
        Vectnfiles = 4;             % Files to open with different levels of colocalization
        VectnChannels=[3:3:21,36];  % Number of channels              
        ChoosingMethod = [1,2,3,5,13]; % Number of methods to compare        
        %%%%%%%%%%%%%%%%%%%%%%%-Parameters-%%%%%%%%%%%%%%%%%%%%%%%
end


% Parameter definiton
theta1 = 0.01;
theta2 = 0.01;
alphaA = 0.001;
alphaX = 0.001;
% Directory to load the H concentration matrix from .mat. With different
% levels of colocalization.
SampleDataFolder = 'Synthetic data/';
% Directory to generate/load the theoretical Spectra. It is used as
% ground_truth
Spectra = 'AlexaRPEPEC'; % Choose which spectra do you want to test. 
TheoreticalSpectraFolder = 'd2_A_from_TheorSpect_';
addpath(genpath('TheoreticalSpectra'));

% Load names of files, and number of files.
imagefiles = dir([SampleDataFolder,'*.mat']);      
nfiles = length(imagefiles);

num = 1;

errorCounting = 1;
ErrortoHistogram=[];
mixHerrorTotal=[];
HmixingErrorAcumlteTotal=[];

RACC_NMFML=[]; Time_NMFML=[];
RACC_NMFML_Exp=[]; Time_NMFML_Exp=[];
RACC_NMFSB=[]; Time_NMFSB=[];
RACC_Names=[]; Time_Names=[];
RACC_Theoretical=[]; Time_Theoretical=[];
RACC_NMFSB_Gauss=[]; Time_NMFSB_Gauss=[];
RACC_NMFRI=[]; Time_NMFRI=[];
RACC_NMFRI_Exp=[]; Time_NMFRI_Exp=[];
RACC_NMFSB_WithoutPixlSel=[]; Time_NMFSB_WithoutPixlSel=[];
RACC_NMFML_WithoutPixlSel=[]; Time_NMFML_WithoutPixlSel=[];
RACC_NMFRI_WithoutPixlSel=[]; Time_NMFRI_WithoutPixlSel=[];
RACC_GaussianInit=[]; Time_GaussianInit=[];
RACC_ExponentialInit=[]; Time_ExponentialInit=[];
RACC_NMFTBCR=[]; Time_NMFTBCR=[];
RACC_NMFTBCR_Exp=[]; Time_NMFTBCR_Exp=[];
RACC_NMFTBCR_WithoutPixlSel=[]; Time_NMFTBCR_WithoutPixlSel=[];
RACC_NMFTBCR_GraphCut = []; Time_NMFTBCR_GraphCut = [];
RACC_IPNMF_VCA = []; Time_IPNMF_VCA = [];
RACC_IPNMF = []; Time_IPNMF = [];


for RMSEThRef=VectRMSEThRef        
for nc=VectnChannels 
    % Creating A mixing matrix with user-defined number of channels.
    Aref = feval(str2func(strcat(TheoreticalSpectraFolder,Spectra)), nc);
    
    % Sum-to-one normalization
    Aref(Aref==0)=100*eps;
    Aref = Aref*diag(1./(sum(Aref,1)+eps));

    % Load Calculated spectra -A0- (using different levels of RMSEthRef), or calculate ir in each iteration.
    % We load them because it has to be the same as the used in the Neher
    % Plugin.
    % Calculate the initialization matrix applying additive and multiplicative
    % factors to the reference spectra.
    LoadSpectraCalculated=true;
    if LoadSpectraCalculated
        if length(VectRMSEThRef)>1 
            if RMSEThRef ==2 
                load('TheoreticalSpectra/A0_changingSimilarity_200.mat')
            elseif RMSEThRef ==2.5 
                load('TheoreticalSpectra/A0_changingSimilarity_250.mat')                  
            elseif RMSEThRef ==3 
                load('TheoreticalSpectra/A0_changingSimilarity_300.mat')            
            elseif RMSEThRef ==3.5 
                load('TheoreticalSpectra/A0_changingSNR.mat')
            elseif RMSEThRef ==4.25 
                load('TheoreticalSpectra/A0_changingSimilarity_425.mat')
            elseif RMSEThRef ==4 
                load('TheoreticalSpectra/A0_changingSimilarity_400.mat')
            elseif RMSEThRef ==4.5 
                load('TheoreticalSpectra/A0_changingSimilarity_450.mat')
            else
                rsum = (rand(size(Aref))-0.5)*(0.0005*(RMSEThRef-1)); % additive factor 
                rmult=(rand(size(Aref))-0.5)*(2.5*(RMSEThRef-1)); rmult(rmult<=0) = 0; % multiplicative factor.   
                A0 = (Aref+ rsum).*(1+rmult);  
                A0=smoothdata(A0,1,'SmoothingFactor',0.1);
                A0(A0<=0)=0.0001; 
                A0 = A0*diag(1./(sum(A0,1)+eps));

%             elseif RMSEThRef ==5 
%                 load('TheoreticalSpectra/A0_changingSimilarity_500.mat')
            end
        end
    else 
        rsum = (rand(size(Aref))-0.5)*(0.0005*(RMSEThRef-1)); % additive factor 
        rmult=(rand(size(Aref))-0.5)*(2.5*(RMSEThRef-1)); rmult(rmult<=0) = 0; % multiplicative factor.   
        A0 = (Aref+ rsum).*(1+rmult);  
        A0=smoothdata(A0,1,'SmoothingFactor',0.1);
        A0(A0<=0)=0.0001; 
        A0 = A0*diag(1./(sum(A0,1)+eps));
    end
    
    % Colocalization study
    if length(Vectnfiles)>1 
        load('A0_changingSNR.mat')
    end
    
    % SNR study
    if length(VectSNR)>1 
        load('A0_changingSNR.mat')
    end
    
    % Spectral Resolution Study. 
    if length(VectnChannels)>1
        % From one channel choose 
        load('A0_changingChannels.mat')                       
        % Create ultra high resolution spectra.
        for A0i=1:size(A0,2)
            newA0(:,A0i) = interp1(1:size(A0,1),A0(:,A0i),0.1:0.1:size(A0,1));
        end        
        newA0(isnan(newA0))=0.0001;        
        % Downsample the ultra high resolution spectra to specified number of channels.
        for newA0i=1:nc
            newnewA0(newA0i,:) =  sum(newA0((newA0i-1)*size(newA0,1)/nc+1:newA0i*size(newA0,1)/nc,:));
        end
        newnewA0(newnewA0<=0)=0.0001; 
        A0 = newnewA0*diag(1./(sum(newnewA0,1)+eps));
        clear newA0;
        clear newnewA0;
    end


for file=Vectnfiles 
    
    %Load concentration matrix, Ground Truth, H.
    currentfilename = imagefiles(file).name;
    GT = load([SampleDataFolder,currentfilename]);
    GT = GT.GT;
    disp(num2str(currentfilename));
    
    % Avoid negative concent Carpetarations due to simulated image noise on
    % H.
    GT(GT<0) = eps;
    GT = GT(:,1:length(GT)/10000:end);
    
for SNR =VectSNR 
      
    % Creation of mixing matrix simulating filtering of Y to create Hfilt
    AFilt = [ones(nc/3,1),zeros(nc/3,1),zeros(nc/3,1); ...
                 zeros(nc/3,1),ones(nc/3,1),zeros(nc/3,1); ...
                 zeros(nc/3,1),zeros(nc/3,1),ones(nc/3,1)];
    AFilt = AFilt*diag(1./(sum(AFilt,1)+eps));
    
    % Calculate Noise, N, to add to Y.
    SigStd = mean(mean(Aref*GT))/(10^(SNR/20));
    NoiseSignal = awgn((Aref)*GT,SNR,'measured');
    NoiseSignal = NoiseSignal;
    DarkCurrent = repmat(abs(awgn(mean((Aref)*GT),SNR,'measured')),[size(Aref,1),1]);
    
    % Background pixels are added.
    Background = ones(size(NoiseSignal))*50;
    Background = awgn(Background,10,'measured'); Background = Background./sum(Background);
    
    % Creation of Y matrix with stained pixels and bakground.
    Y = [NoiseSignal + DarkCurrent,Background + DarkCurrent,Background + DarkCurrent];
    Y = Y(:,1:1:end);
    
    % Data Preprocessing.
    [Y, Y_sps, H_sps] = dataPreprocessing(A0,Y);
    
    % Simplex projection of data.
    [UU, SS, VV] = svds(Y./sum(Y),size(A0,2));
    LowHaux = UU'*(Y./sum(Y));
    figure; plot(LowHaux(2,:),LowHaux(3,:),'b.'); hold on;%ALL pixels
    %LowPure = UU'*Y_purepixels; plot(LowPure(2,:),LowPure(3,:),'c.'); %Pure Pixels
    LowSparsest = UU'*(Y_sps./sum(Y_sps)); plot(LowSparsest(2,:),LowSparsest(3,:),'c.'); %Sparsest Pixels
    LowAnorm = UU'*Aref; hold on; plot(LowAnorm(2,:),LowAnorm(3,:),'go','MarkerSize',7); %Reference Spectra
    LowA0 = UU'*A0; hold on; plot(LowA0(2,:),LowA0(3,:),'k+','MarkerSize',7); %Initial Spectra
    xlabel('Subspace Vector 3'); ylabel('Subspace Vector 1'); title('Simplex Projection')        

for method = ChoosingMethod %To calculate Anmf.
        switch method 
            case 1 % NMF-ML using theoretical Spectra                               
                
                [A, H,~,~, elapsedTime] = NMF_ML(Y_sps,[],'A0',A0,'H0',H_sps);
            case 2 % NMF-SB using theoretical Spectra
                % Save Observation matrix, so that it can be opened using NeherPlugin. 
                YNeher = uint16(Y_sps./max(Y_sps(:)).*65535);
                imwrite(YNeher(1,:),['YtoNeher3/Y_AlexaRPEPEC_Col',num2str(file),'_Similarity',num2str(RMSEThRef*100),'_SNR',num2str(SNR),'_nChnnls',num2str(nc),'.tif'])
                for i =2:size(YNeher,1)
                    imwrite(YNeher(i,:),['YtoNeher3/Y_AlexaRPEPEC_Col',num2str(file),'_Similarity',num2str(RMSEThRef*100),'_SNR',num2str(SNR),'_nChnnls',num2str(nc),'.tif'],'WriteMode','append')
                end                                                                                                                                                 

                clear Aneher;
                % Load calculated spectra from folder
                fnameS = ['YtoNeher3/Y_AlexaRPEPEC_Col',num2str(file),'_Similarity',num2str(RMSEThRef*100),'_SNR',num2str(SNR),'_nChnnls',num2str(nc),'_Spectrum1.txt'];
                fSpectraNEHER = load(fnameS); Aneher(:,1) = fSpectraNEHER(:,2);
                fnameS = ['YtoNeher3/Y_AlexaRPEPEC_Col',num2str(file),'_Similarity',num2str(RMSEThRef*100),'_SNR',num2str(SNR),'_nChnnls',num2str(nc),'_Spectrum2.txt'];
                fSpectraNEHER = load(fnameS); Aneher(:,2) = fSpectraNEHER(:,2);
                fnameS = ['YtoNeher3/Y_AlexaRPEPEC_Col',num2str(file),'_Similarity',num2str(RMSEThRef*100),'_SNR',num2str(SNR),'_nChnnls',num2str(nc),'_Spectrum3.txt'];
                fSpectraNEHER = load(fnameS); Aneher(:,3) = fSpectraNEHER(:,2);
                A = Aneher*diag(1./(sum(Aneher,1)+eps));
                elapsedTime = 60;

            case 3 % Theoretical Spectra                
                A = A0;           
            case 4 % NMF-SB using gaussian initialization.
                % Do the Neher Thing
                YNeher = uint16(Y_sps./max(Y_sps(:)).*65535);
                imwrite(YNeher(1,:),['YtoNeher3Gauss/Y_AlexaRPEPEC_Col',num2str(file),'_Similarity',num2str(RMSEThRef*100),'_SNR',num2str(SNR),'_nChnnls',num2str(nc),'.tif'])
                for i =2:size(YNeher,1)
                    imwrite(YNeher(i,:),['YtoNeher3Gauss/Y_AlexaRPEPEC_Col',num2str(file),'_Similarity',num2str(RMSEThRef*100),'_SNR',num2str(SNR),'_nChnnls',num2str(nc),'.tif'],'WriteMode','append')
                end                                                                                                                        

                clear Aneher;
                % Load calculated spectra from folder
                fnameS = ['YtoNeher3Gauss/Y_AlexaRPEPEC_Col',num2str(file),'_Similarity',num2str(RMSEThRef*100),'_SNR',num2str(SNR),'_nChnnls',num2str(nc),'_spectrum1.txt'];
                fSpectraNEHER = load(fnameS); Aneher(:,1) = fSpectraNEHER(:,2);
                fnameS = ['YtoNeher3Gauss/Y_AlexaRPEPEC_Col',num2str(file),'_Similarity',num2str(RMSEThRef*100),'_SNR',num2str(SNR),'_nChnnls',num2str(nc),'_spectrum2.txt'];
                fSpectraNEHER = load(fnameS); Aneher(:,2) = fSpectraNEHER(:,2);
                fnameS = ['YtoNeher3Gauss/Y_AlexaRPEPEC_Col',num2str(file),'_Similarity',num2str(RMSEThRef*100),'_SNR',num2str(SNR),'_nChnnls',num2str(nc),'_spectrum3.txt'];
                fSpectraNEHER = load(fnameS); Aneher(:,3) = fSpectraNEHER(:,2);
                A = Aneher*diag(1./(sum(Aneher,1)+eps));
                elapsedTime =60;
            
            case 5 % NMF-RI using theoretical Spectra                                  
                [A, H, elapsedTime] =  NMF_RI(Y_sps,A0,H_sps, theta1, theta2, alphaA, alphaX);
            
            case 6 % NMF-ML using Exponential Spectra
                % In this case we recalculate the H
                H0exp = max(100*eps,pinv(exponential_matrix(nf,nc))*Y_sps); 
                [A, H,~,~,elapsedTime] = NMF_ML(Y_sps,[],'A0',exponential_matrix(nf,nc),'H0',H0exp);
            
            case 7 % NMF-RI using Exponential Spectra
                H0gauss = max(100*eps,pinv(exponential_matrix(nf,nc))*Y_sps); 
                [A, H,elapsedTime] =  NMF_RI(Y_sps,exponential_matrix(nf,nc),H0gauss, theta1, theta2, alphaA, alphaX);
            
            case 8 % Neher using Theoretical Spectra without Sparseness pixelSelection
                % Do the Neher Thing
                YNeher = uint16(Y./max(Y(:)).*65535);
                imwrite(YNeher(1,:),['YtoNeher3NoPixlSel/Y_AlexaRPEPEC_Col',num2str(file),'_Similarity',num2str(RMSEThRef*100),'_SNR',num2str(SNR),'_nChnnls',num2str(nc),'.tif'])
                for i =2:size(YNeher,1)
                    imwrite(YNeher(i,:),['YtoNeher3NoPixlSel/Y_AlexaRPEPEC_Col',num2str(file),'_Similarity',num2str(RMSEThRef*100),'_SNR',num2str(SNR),'_nChnnls',num2str(nc),'.tif'],'WriteMode','append')
                end                                                                                                                                                 

                clear Aneher;
                % Load calculated spectra from folder
                fnameS = ['YtoNeher3NoPixlSel/Y_AlexaRPEPEC_Col',num2str(file),'_Similarity',num2str(RMSEThRef*100),'_SNR',num2str(SNR),'_nChnnls',num2str(nc),'_Spectrum1.txt'];
                fSpectraNEHER = load(fnameS); Aneher(:,1) = fSpectraNEHER(:,2);
                fnameS = ['YtoNeher3NoPixlSel/Y_AlexaRPEPEC_Col',num2str(file),'_Similarity',num2str(RMSEThRef*100),'_SNR',num2str(SNR),'_nChnnls',num2str(nc),'_Spectrum2.txt'];
                fSpectraNEHER = load(fnameS); Aneher(:,2) = fSpectraNEHER(:,2);
                fnameS = ['YtoNeher3NoPixlSel/Y_AlexaRPEPEC_Col',num2str(file),'_Similarity',num2str(RMSEThRef*100),'_SNR',num2str(SNR),'_nChnnls',num2str(nc),'_Spectrum3.txt'];
                fSpectraNEHER = load(fnameS); Aneher(:,3) = fSpectraNEHER(:,2);
                A = Aneher*diag(1./(sum(Aneher,1)+eps));
                elapsedTime = 60;

            case 9 % NMF-ML without PixelSelection method
                [A, H,~,~,elapsedTime] = NMF_ML(Y,[],'A0',A0,'H0',H0);

            case 10 % NMF-RI without PixelSelection method
                [A, H,elapsedTime] =  NMF_RI(Y,A0,H0, theta1, theta2, alphaA, alphaX);
            
            case 11 % Gaussian Initialization
                A = gauss_matrix(nf,nc);
            
            case 12 % Exponential Initialization
                A = exponential_matrix(nf,nc);
                
            case 13 % TBCR-NMF using theoretical spectra
                addpath('TBCR-NMF')
                phi = 0.1; 
                [C,S,time,cost_record,lambda_sq,elapsedTime]= mynmf_ghals_L12_Ultimate(Y_sps',H_sps',A0',phi); 
                A = S'; 
            case 14 % TBCR-NMF using Exponential Spectra
                addpath('TBCR-NMF')
                phi = 0.1; 
                H0exp = max(100*eps,pinv(exponential_matrix(nf,nc))*Y_sps);                
                [C,S,time,cost_record,lambda_sq, elapsedTime]= mynmf_ghals_L12_Ultimate(Y_sps',H0exp',exponential_matrix(nf,nc)',phi);
                A = S'; 
            case 15 % Without pixel Selection 
                addpath('TBCR-NMF')
                phi = 0.1; 
                [C,S,time,cost_record,lambda_sq, elapsedTime]= mynmf_ghals_L12_Ultimate(Y',H0',A0',phi);
                A = S'; 
            case 16 % TBCR with Graph-Cut pixel selection
                addpath('TBCR-NMF')
                phi = 0.1; 
                ACutInit = graphCutInit(Y_sps,size(A0,2));                
                [C,S,time,cost_record,lambda_sq, elapsedTime]= mynmf_ghals_L12_Ultimate(Y_sps',H_sps',ACutInit',phi);
                A = S'; 
            case 17 % IP-NMF with theoretical spectra
                addpath('IP-NMF')
                % Initialization
                M = size(A0,2);
                L = size(A0,1);
                X = Y'; 
                P = size(X,1);
                R_ref_VCA = NaN(M,L);
                R_ref_VCA(:,~isnan(X(1,:))) = A0';            
                R_ref_VCA = repmat(R_ref_VCA,P,1);
                R0 = R_ref_VCA;     
                
                % C
                C0_IPNMF = zeros(P,P*M);
                for p = 1 : P
                    C0_IPNMF(p,(p-1)*M+1:p*M) = FCLSU(X(p,:)',R0(p,:)');
                end                  
                [Cp_t,Rp_t,Iter,Cp1,Rp1] = Grad_proj_pixels_contraintes_V2(X, C0_IPNMF, R_ref_VCA,3);
                for p = 1 : P
                    R_p = Rp_t((p-1)*M+1:p*M,:);
                    Hnmf(:,p) = max(100*eps,pinv(R_p')*Y(:,p)); 
                end                                 
                               
            case 18 % IP-NMF with VCA initialization
                addpath('IP-NMF')
%                 Y_sps = Y_sps(:,1:100:end);
                
                % Initialization
                Results_VCA = EIA_VCA(Y,1,size(A0,2),false);            
                M = size(A0,2);
                L = size(A0,1);
                X = Y'; 
                P = size(X,1);
                R_ref_VCA = NaN(M,L);
                R_ref_VCA(:,~isnan(X(1,:))) = Results_VCA.E';            
                R_ref_VCA = repmat(R_ref_VCA,P,1);
                R0 = R_ref_VCA;                

                C0_IPNMF = zeros(P,P*M);
                for p = 1 : P
                    C0_IPNMF(p,(p-1)*M+1:p*M) = FCLSU(X(p,:)',R0(p,:)');
                end                                   
                [Cp_t,Rp_t,Iter,C1,R1] = Grad_proj_pixels_contraintes_V2(X, C0_IPNMF, R_ref_VCA,3);
                
                for p = 1 : P
                    R_p = Rp_t((p-1)*M+1:p*M,:);
                    Hnmf(:,p) = max(100*eps,pinv(R_p')*Y(:,p)); 
                end                                 
            end
        
        % Sum-to-one normalization of every Blind spectral unmixing method.
        A = A*diag(1./(sum(A,1)+eps));

        if method<17 % For IPNMF method use the alternative calculation of the Hnmf
            % Calculate the resulted Hnmf, concentration matrix.
            Hnmf = max(100*eps,pinv(A)*Y); 
        end
        % Calculate the resulted Href, concentration matrix.
        Href = max(100*eps,pinv(Aref)*Y);
        
        if length(VectnChannels)>1  
            Href = [GT,Hnmf(:,size(GT,2)+1:end)];
        end
        
        % Concentration matrix using the initialization spectra
        H0 = max(100*eps,pinv(A0)*Y);
                        
        % Creating concentration matrix, using the filter method.
        Hfilt = max(100*eps,pinv(AFilt)*Y);
        
        % Used to calculate the RMSE between the theoretical and the
        % reference spectra 
        if method==1
            if ~exist('RMSEworst','var') 
                RMSEworst = sqrt(mean((Hnmf(:)-H0(:)).^2));
            else
                RMSEworst = [RMSEworst, sqrt(mean((Hnmf(:)-H0(:)).^2))];
            end                        
        end
        A0 = A0*diag(1./(sum(A0,1)+eps));
        
        % This is used to print the extracted spectra.
        LowA = UU'*A; hold on; plot(LowA(2,:),LowA(3,:),'r+','MarkerSize',7);
        legend('Non segmented Pixels','Segmented Pixels','Reference Spectra','Initialization Spectra','Calculated Spectra')

        % RMSE between Href and Hnmf
        RMSErefNMF(errorCounting) = sqrt(mean((Href(:)-Hnmf(:)).^2));                    
        
        % RMSE between Href and Hfilt
        RMSErefFilt(errorCounting) = sqrt(mean((Href(:)-Hfilt(:)).^2));                        
        
        % Calculate RACC, Performance evaluation metric.
        if RMSErefNMF(errorCounting)<RMSErefFilt(errorCounting)             
            RACC(errorCounting) = (RMSErefFilt(errorCounting)-RMSErefNMF(errorCounting))/RMSErefFilt(errorCounting);
        else            
            RACC(errorCounting) = (RMSErefFilt(errorCounting)-RMSErefNMF(errorCounting))/RMSErefNMF(errorCounting);
        end        
        disp(['RelativeAccuracy using ',num2str(nc),' Channels: ',num2str(RACC(errorCounting) )]);
        RACC_Names = [RACC_Names,{[currentfilename,'_Chnnl',num2str(nc),'_Spectra',Spectra]}];
        % Values that are negative put to zero.
        RACC(RACC<0)=0;
        switch method            
            case 1
                 RACC_NMFML = [RACC_NMFML, RACC(errorCounting)];                 
                 Time_NMFML = [Time_NMFML, elapsedTime];                 
            case 2
                 RACC_NMFSB = [RACC_NMFSB, RACC(errorCounting)];   
                 Time_NMFSB = [Time_NMFSB, elapsedTime];   
            case 3
                 RACC_Theoretical = [RACC_Theoretical, RACC(errorCounting)];   
                 Time_Theoretical = [Time_Theoretical, elapsedTime];   
            case 4
                 RACC_NMFSB_Gauss = [RACC_NMFSB_Gauss, RACC(errorCounting)];   
                 Time_NMFSB_Gauss = [Time_NMFSB_Gauss, elapsedTime];   
            case 5
                 RACC_NMFRI = [RACC_NMFRI, RACC(errorCounting)];   
                 Time_Theoretical = [Time_Theoretical, elapsedTime];   
            case 6 
                 RACC_NMFML_Exp = [RACC_NMFML_Exp, RACC(errorCounting)];   
                 Time_NMFML_Exp = [Time_NMFML_Exp, elapsedTime];   
            case 7 
                RACC_NMFRI_Exp = [RACC_NMFRI_Exp, RACC(errorCounting)];  
                Time_NMFRI_Exp = [Time_NMFRI_Exp, elapsedTime];  
            case 8
                RACC_NMFSB_WithoutPixlSel = [RACC_NMFSB_WithoutPixlSel, RACC(errorCounting)];  
                Time_NMFRI_Exp = [Time_NMFRI_Exp, elapsedTime];  
            case 9
                RACC_NMFML_WithoutPixlSel = [RACC_NMFML_WithoutPixlSel, RACC(errorCounting)];  
                Time_NMFRI_Exp = [Time_NMFRI_Exp, elapsedTime];  
            case 10
                RACC_NMFRI_WithoutPixlSel = [RACC_NMFRI_WithoutPixlSel, RACC(errorCounting)];
                Time_NMFRI_WithoutPixlSel = [Time_NMFRI_WithoutPixlSel, elapsedTime];
            case 11
                RACC_GaussianInit = [RACC_GaussianInit, RACC(errorCounting)];
                Time_GaussianInit = [Time_GaussianInit, elapsedTime];
            case 12
                RACC_ExponentialInit = [RACC_ExponentialInit, RACC(errorCounting)];
                Time_ExponentialInit = [Time_ExponentialInit, elapsedTime];
            case 13 
                RACC_NMFTBCR = [RACC_NMFTBCR, RACC(errorCounting)]; 
                Time_NMFTBCR = [Time_NMFTBCR, elapsedTime]; 
            case 14
                RACC_NMFTBCR_Exp = [RACC_NMFTBCR_Exp, RACC(errorCounting)]; 
                Time_NMFTBCR_Exp = [Time_NMFTBCR_Exp, elapsedTime]; 
            case 15
                RACC_NMFTBCR_WithoutPixlSel = [RACC_NMFTBCR_WithoutPixlSel, RACC(errorCounting)]; 
                Time_NMFTBCR_WithoutPixlSel = [Time_NMFTBCR_WithoutPixlSel, elapsedTime]; 
            case 16
                RACC_NMFTBCR_GraphCut = [RACC_NMFTBCR_GraphCut, RACC(errorCounting)]; 
                Time_NMFTBCR_GraphCut = [Time_NMFTBCR_GraphCut, elapsedTime]; 
            case 17
                RACC_IPNMF= [RACC_IPNMF, RACC(errorCounting)]; 
                Time_IPNMF= [Time_IPNMF, elapsedTime]; 
            case 18
                RACC_IPNMF_VCA = [RACC_IPNMF_VCA, RACC(errorCounting)]; 
                Time_IPNMF_VCA = [Time_IPNMF_VCA, elapsedTime]; 
        end        
        errorCounting = errorCounting+1;        
end

% Accumulate in one matrix the results of the methods used. 
errorCounting = 1;
mixHerror(method,:) = RMSErefNMF;
HmixingErrorAcumlte(method,:) = RMSErefFilt;
RelativeAccuracyAcumlte(method,:) = RACC;
end
    % Accumulate in one matrix the results of the methods used. 
    if isempty(mixHerrorTotal) %Initialize mixHerrorTotal
        mixHerrorTotal = mixHerror;
        HmixingErrorAcumlteTotal = HmixingErrorAcumlte;
        RelativeAccuracyAcumlteTotal = RelativeAccuracyAcumlte;
    else
        mixHerrorTotal = mixHerrorTotal + mixHerror;
        HmixingErrorAcumlteTotal = HmixingErrorAcumlteTotal + HmixingErrorAcumlte;
        RelativeAccuracyAcumlteTotal = RelativeAccuracyAcumlteTotal + RelativeAccuracyAcumlte; 
    end

end
end
end

%% Display Results
% Show A matrices
figure;
vectorA = 1:size(Aref,1);
plot(vectorA,Aref(:,1),'Color','g','LineStyle','-','LineWidth',2);
hold on;
plot(vectorA,Aref(:,2),'Color','b','LineStyle','-','LineWidth',2);
plot(vectorA,Aref(:,3),'Color','r','LineStyle','-','LineWidth',2);
plot(vectorA,A(:,1),'Color','g','Marker','o','LineStyle','-','LineWidth',1);
plot(vectorA,A(:,2),'Color','b','Marker','o','LineStyle','-','LineWidth',1);
plot(vectorA,A(:,3),'Color','r','Marker','o','LineStyle','-','LineWidth',1);
plot(vectorA,A0(:,1),'Color','g','LineStyle','--','LineWidth',2);
plot(vectorA,A0(:,2),'Color','b','LineStyle','--','LineWidth',2);
plot(vectorA,A0(:,3),'Color','r','LineStyle','--','LineWidth',2);
legend('Alexa488 Reference','R-PE Reference','PE-cy5 Reference','Alexa488 NMF-RI','R-PE NMF-RI','PE-cy5 NMF-RI','Alexa488 Initialization','R-PE Initialization','PE-cy5 Initialization');
xlabel('Spectral Channels')
ylabel('Normalized intensity')
axis tight
title('Comparing mixing matrices')


switch Study
    case 'Colocalization'
        % Display Spatial Colocalization
        ColLevels = [0,25,50,75,100];
        figure; 
        plot(ColLevels,ones(size(ColLevels)),'LineWidth',2,'LineStyle','-');hold on;
        plot(ColLevels,zeros(size(ColLevels)),'LineWidth',2,'LineStyle','-');hold on;
        plot(ColLevels,RACC_Theoretical,'LineWidth',2,'LineStyle','-'); 
        plot(ColLevels,RACC_NMFML,'LineWidth',2,'LineStyle','-'); 
        plot(ColLevels,RACC_NMFSB,'LineWidth',2,'LineStyle','-'); 
        plot(ColLevels,RACC_NMFRI,'LineWidth',2,'LineStyle','-');
        plot(ColLevels,RACC_NMFTBCR,'LineWidth',2,'LineStyle','-');
%         plot(ColLevels,RACC_IPNMF,'LineWidth',2,'LineStyle','-','Color',[0.25, 0.25, 0.25]);        
        plot(ColLevels,RACC_NMFML_Exp,'LineWidth',2,'LineStyle','--','Color',[0.4940, 0.1840, 0.5560]);
        plot(ColLevels,RACC_NMFSB_Gauss,'LineWidth',2,'LineStyle','--','Color',[0.4660, 0.6740, 0.1880]);
%         plot(ColLevels,RACC_NMFRI_Exp,'LineWidth',2,'LineStyle','--','Color',[0.3010, 0.7450, 0.9330]);
        plot(ColLevels,RACC_NMFTBCR_GraphCut,'LineWidth',2,'LineStyle','--','Color',[0.6350, 0.0780, 0.1840]);        
%         plot(ColLevels,RACC_IPNMF_VCA,'LineWidth',2,'LineStyle','--','Color',[0.25, 0.25, 0.25]);        
        % plot(ColLevels,RelAcctoHistogram_Neher_WithoutPixlSel,'LineWidth',2,'LineStyle','--','Color',[0.4660, 0.6740, 0.1880]);
        % plot(ColLevels,RelAcctoHistogram_ML_WithoutPixlSel,'LineWidth',2,'LineStyle','--','Color',[0.4940, 0.1840, 0.5560]);
        % plot(ColLevels,RelAcctoHistogram_RI_WithoutPixlSel,'LineWidth',2,'LineStyle','--','Color',[0.3010, 0.7450, 0.9330]);
%         legend('Reference Spectra','Filter method','Theoretical Spectra','NMF-ML','NMF-SB','NMF-RI','TBCR-NMF','IP-NMF','NMF-ML(Exponential Initialization)','NMF-SB(Gaussian Initialization)','NMF-RI(Exponential Initialization)','TBCR-NMF(nCut Initialization)','IP-NMF(VCA Initialization)');
        legend('Reference Spectra','Filter method','Theoretical Spectra','NMF-ML','NMF-SB','NMF-RI','TBCR-NMF','NMF-ML(Exponential Initialization)','NMF-SB(Gaussian Initialization)','TBCR-NMF(nCut Initialization)');
        xlabel('Degrees of spatial colocalization (%)');
        ylabel('RACC');
        xticks([0,25,50,75,100]); xticklabels({'0','25','50','75','100'});
        
    case 'RMSEThRef'
        % Display RMSE between Theoretical and Reference Spectra
        [RMSEworstSort,RMSEworstSortIndex] = sort(RMSEworst); % sort element from less to more RMSE.
        figure; 
        plot(RMSEworstSort,ones(size(RACC_NMFML(RMSEworstSortIndex))),'LineWidth',2,'LineStyle','-');hold on;
        plot(RMSEworstSort,zeros(size(RACC_NMFML(RMSEworstSortIndex))),'LineWidth',2,'LineStyle','-');hold on;
        plot(RMSEworstSort,RACC_Theoretical(RMSEworstSortIndex),'LineWidth',2,'LineStyle','-'); 
        plot(RMSEworstSort,RACC_NMFML(RMSEworstSortIndex),'LineWidth',2,'LineStyle','-','Color',[0.4940, 0.1840, 0.5560]); 
        plot(RMSEworstSort,RACC_NMFSB(RMSEworstSortIndex),'LineWidth',2,'LineStyle','-','Color',[0.4660, 0.6740, 0.1880]); 
        plot(RMSEworstSort,RACC_NMFRI(RMSEworstSortIndex),'LineWidth',2,'LineStyle','-','Color',[0.3010, 0.7450, 0.9330]);
        plot(RMSEworstSort,RACC_NMFTBCR(RMSEworstSortIndex),'LineWidth',2,'LineStyle','-','Color',[0.6350, 0.0780, 0.1840]);%         
        legend('Reference Spectra','Filter method','Theoretical Spectra','NMF-ML','NMF-SB','NMF-RI','TBCR-NMF');
        xlabel('RMSE between Theoretical and Reference Spectra');
        ylabel('RACC'); axis tight
%         xticks(RMSEworstSort); xticklabels({'0','25','50','75','100'});

    case 'SNR'
        iter = 1;
        RACC_Theoretical_reordered=zeros(length(VectSNR),1);RACC_NMFML_reordered=zeros(length(VectSNR),1);
        RACC_NMFSB_reordered=zeros(length(VectSNR),1); RACC_NMFRI_reordered=zeros(length(VectSNR),1);
        RACC_NMFTBCR_reordered=zeros(length(VectSNR),1);RACC_NMFSB_WithoutPixlSel_reordered=zeros(length(VectSNR),1);
        RACC_NMFML_WithoutPixlSel_reordered=zeros(length(VectSNR),1); RACC_NMFRI_WithoutPixlSel_reordered=zeros(length(VectSNR),1);
        RACC_NMFTBCR_WithoutPixlSel_reordered=zeros(length(VectSNR),1);
        for file=1:length(Vectnfiles)
            for SNR =1:length(VectSNR) 
                RACC_Theoretical_reordered(SNR) = RACC_Theoretical_reordered(SNR) + RACC_Theoretical(iter);
                RACC_NMFML_reordered(SNR) = RACC_NMFML_reordered(SNR) + RACC_NMFML(iter);
                RACC_NMFSB_reordered(SNR) = RACC_NMFSB_reordered(SNR) + RACC_NMFSB(iter);
                RACC_NMFRI_reordered(SNR) = RACC_NMFRI_reordered(SNR) + RACC_NMFRI(iter);
                RACC_NMFTBCR_reordered(SNR) = RACC_NMFTBCR_reordered(SNR) + RACC_NMFTBCR(iter);
                RACC_NMFSB_WithoutPixlSel_reordered(SNR) = RACC_NMFSB_WithoutPixlSel_reordered(SNR) + RACC_NMFSB_WithoutPixlSel(iter);
                RACC_NMFTBCR_WithoutPixlSel_reordered(SNR) = RACC_NMFTBCR_WithoutPixlSel_reordered(SNR) + RACC_NMFTBCR_WithoutPixlSel(iter);
                RACC_NMFRI_WithoutPixlSel_reordered(SNR) = RACC_NMFRI_WithoutPixlSel_reordered(SNR) + RACC_NMFRI_WithoutPixlSel(iter);
                RACC_NMFML_WithoutPixlSel_reordered(SNR) = RACC_NMFML_WithoutPixlSel_reordered(SNR) + RACC_NMFML_WithoutPixlSel(iter);
                iter = iter+1;
            end
        end
        RACC_Theoretical_reordered=RACC_Theoretical_reordered/length(Vectnfiles);
        RACC_NMFML_reordered = RACC_NMFML_reordered/length(Vectnfiles);
        RACC_NMFSB_reordered = RACC_NMFSB_reordered/length(Vectnfiles);
        RACC_NMFRI_reordered = RACC_NMFRI_reordered/length(Vectnfiles);
        RACC_NMFTBCR_reordered = RACC_NMFTBCR_reordered/length(Vectnfiles);
        RACC_NMFSB_WithoutPixlSel_reordered = RACC_NMFSB_WithoutPixlSel_reordered/length(Vectnfiles);
        RACC_NMFTBCR_WithoutPixlSel_reordered = RACC_NMFTBCR_WithoutPixlSel_reordered/length(Vectnfiles);
        RACC_NMFRI_WithoutPixlSel_reordered = RACC_NMFRI_WithoutPixlSel_reordered/length(Vectnfiles);
        RACC_NMFML_WithoutPixlSel_reordered = RACC_NMFML_WithoutPixlSel_reordered/length(Vectnfiles);
        % SNR 
        figure; 
        plot(VectSNR,ones(size(VectSNR)),'LineWidth',2,'LineStyle','-');hold on;
        plot(VectSNR,zeros(size(VectSNR)),'LineWidth',2,'LineStyle','-');hold on;
        plot(VectSNR,RACC_Theoretical_reordered,'LineWidth',2,'LineStyle','-'); 
        plot(VectSNR,RACC_NMFML_reordered,'LineWidth',2,'LineStyle','-'); 
        plot(VectSNR,RACC_NMFSB_reordered,'LineWidth',2,'LineStyle','-'); 
        plot(VectSNR,RACC_NMFRI_reordered,'LineWidth',2,'LineStyle','-');
        plot(VectSNR,RACC_NMFTBCR_reordered,'LineWidth',2,'LineStyle','-');
%         plot(VectSNR,RACC_IPNMF,'LineWidth',2,'LineStyle','-','Color', [0.25,0.25,0.25]);
        plot(VectSNR,RACC_NMFSB_WithoutPixlSel_reordered,'LineWidth',2,'LineStyle','--','Color',[0.4660, 0.6740, 0.1880]);
        plot(VectSNR,RACC_NMFML_WithoutPixlSel_reordered,'LineWidth',2,'LineStyle','--','Color',[0.4940, 0.1840, 0.5560]);
        plot(VectSNR,RACC_NMFRI_WithoutPixlSel_reordered,'LineWidth',2,'LineStyle','--','Color',[0.3010, 0.7450, 0.9330]);
        plot(VectSNR,RACC_NMFTBCR_WithoutPixlSel_reordered,'LineWidth',2,'LineStyle','--','Color',[0.6350, 0.0780, 0.1840]);
        legend('Reference Spectra','Filter method','Theoretical Spectra','NMF-ML','NMF-SB','NMF-RI','TBCR-NMF','NMF-SB(Without Pixel Selection)','NMF-ML(Without Pixel Selection)','NMF-RI(Without Pixel Selection)','TBCR-NMF(Without Pixel Selection)');
        xlabel('SNR (dB)');
        ylabel('RACC');
        title('');
        xticks(VectSNR); xticklabels({'20','25','30','35','40'});

    case 'SpectralResolution'
        % Spectral Resolution Study 
        figure; 
        plot(VectnChannels/3,ones(size(VectnChannels)),'LineWidth',2,'LineStyle','-');hold on;
        plot(VectnChannels/3,zeros(size(VectnChannels)),'LineWidth',2,'LineStyle','-');hold on;
        plot(VectnChannels/3,RACC_Theoretical,'LineWidth',2,'LineStyle','-'); 
        plot(VectnChannels/3,RACC_NMFML,'LineWidth',2,'LineStyle','-'); 
        plot(VectnChannels/3,RACC_NMFSB,'LineWidth',2,'LineStyle','-'); 
        plot(VectnChannels/3,RACC_NMFRI,'LineWidth',2,'LineStyle','-');
        plot(VectnChannels/3,RACC_NMFTBCR,'LineWidth',2,'LineStyle','-');
        plot(VectnChannels/3,RACC_IPNMF,'LineWidth',2,'LineStyle','-','Color',[0.25, 0.25, 0.25]);

        legend('Reference Spectra','Filter method','Theoretical Spectra','NMF-ML','NMF-SB','NMF-TI','TBCR-NMF','IP-NMF');
        axis tight;
        xlabel('Ratio Channels/Fluorochromes');
        ylabel('RACC');
        title('');
end  