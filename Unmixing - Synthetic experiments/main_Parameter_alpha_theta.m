% Synthetic test and performance comparison between algorithms (NMF-RI)
% It follows NMF-RI as described in the paper "NMF-RI: Blind spectral unmixing of highly mixed multispectral flow and image cytometry"
% We used emission and excitation data to create theoretical spectra.
% Using synthetic datasets to 
% Written by Daniel Jimenez-Sanchez.

clear all;
close all; 
clc;

% To execute different studies:
    % Choose between -> 'Colocalization', 'RMSEThRef', 'SNR',
    % 'SpectralResolution', 'theta', 'alpha'
    Study = 'theta'; 
    
switch Study
    case 'theta'
        %%%%%%%%%%%%%%%%%%%%%%%-Parameters-%%%%%%%%%%%%%%%%%%%%%%%
        % These parameters studies different levels of colocalization [0,25,50,75 and 100].
        nf = 3;                           % number of fluorochromes        
        VectSNR = 35;                     % SNR levels
        VectRMSEThRef = 3.5;              % Similarity between Reference and theoretical spectra
        Vectnfiles = 4:5;                 % Files to open with different levels of colocalization
        VectnChannels=36;                 % Number of channels              
        Vecttheta1 = [0.5, 0.1, 0.05, 0.01, 0.005, 0.001]; % Stopping Criteria intraLayer
        Vecttheta2 = [0.5, 0.1, 0.05, 0.01, 0.005, 0.001]; % Stopping Criteria interLayer
        VectalphaA = 0.001; 
        VectalphaX = 0.001;
        ChoosingMethod = [5]; % Number of methods to compare       
    case 'alpha'
        %%%%%%%%%%%%%%%%%%%%%%%-Parameters-%%%%%%%%%%%%%%%%%%%%%%%
        % These parameters studies different levels of colocalization [0,25,50,75 and 100].
        nf = 3;                           % number of fluorochromes        
        VectSNR = 35;                     % SNR levels
        VectRMSEThRef = 3.5;              % Similarity between Reference and theoretical spectra
        Vectnfiles = 4:5;                 % Files to open with different levels of colocalization
        VectnChannels=36;                 % Number of channels              
        Vecttheta1 = [0.01]; % Stopping Criteria intraLayer
        Vecttheta2 = [0.01]; % Stopping Criteria interLayer
        VectalphaA = [0.5, 0.25, 0.1, 0.01, 0.001, 0.0001]; 
        VectalphaX = [0.5, 0.25, 0.1, 0.01, 0.001, 0.0001];
        ChoosingMethod = [5]; % Number of methods to compare       
    case 'Colocalization'
        %%%%%%%%%%%%%%%%%%%%%%%-Parameters-%%%%%%%%%%%%%%%%%%%%%%%
        % These parameters studies different levels of colocalization [0,25,50,75 and 100].
        nf = 3;                           % number of fluorochromes        
        VectSNR = 35;                     % SNR levels
        VectRMSEThRef = 3.5;              % Similarity between Reference and theoretical spectra
        Vectnfiles = 1:5;                 % Files to open with different levels of colocalization
        VectnChannels=36;                 % Number of channels              
        ChoosingMethod = [1,2,3,4,5,6,7,13,14,16]; % Number of methods to compare        
        %%%%%%%%%%%%%%%%%%%%%%%-Parameters-%%%%%%%%%%%%%%%%%%%%%%%
    case 'RMSEThRef'
        %%%%%%%%%%%%%%%%%%%%%%%-Parameters-%%%%%%%%%%%%%%%%%%%%%%%
        % These parameters studies different RMSE between Theoretical and Reference Spectra.
        nf = 3;                     % number of fluorochromes        
        VectSNR = 35;               % SNR levels
        VectRMSEThRef = 2:0.5:4;    % Similarity between Reference and theoretical spectra
        Vectnfiles = 4;             % Files to open with different levels of colocalization
        VectnChannels=36;           % Number of channels              
        ChoosingMethod = [1,2,3,5,13]; % Number of methods to compare        
        %%%%%%%%%%%%%%%%%%%%%%%-Parameters-%%%%%%%%%%%%%%%%%%%%%%%
    case 'SNR'
        %%%%%%%%%%%%%%%%%%%%%%%-Parameters-%%%%%%%%%%%%%%%%%%%%%%%
        % These parameters studies different levels of SNR.
        nf = 3;                     % number of fluorochromes        
        VectSNR = [20:5:40];        % SNR levels
        VectRMSEThRef = 3.5;        % Similarity between Reference and theoretical spectra
        Vectnfiles = 4;             % Files to open with different levels of colocalization
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

RACC_NMFML=[];
RACC_NMFML_Exp=[];
RACC_NMFSB=[];
RACC_Names=[];
RACC_Theoretical=[];
RACC_NMFSB_Gauss=[];
RACC_NMFRI=[];iterations_NMFRI=[];
RACC_NMFRI_Exp=[];
RACC_NMFSB_WithoutPixlSel=[];
RACC_NMFML_WithoutPixlSel=[];
RACC_NMFRI_WithoutPixlSel=[];
RACC_GaussianInit=[];
RACC_ExponentialInit=[];
RACC_NMFTBCR=[];
RACC_NMFTBCR_Exp=[];
RACC_NMFTBCR_WithoutPixlSel=[];
RACC_NMFTBCR_GraphCut = [];

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
            end
            if RMSEThRef ==2.5 
                    load('TheoreticalSpectra/A0_changingSimilarity_250.mat')
            end        
            if RMSEThRef ==3 
                    load('TheoreticalSpectra/A0_changingSimilarity_300.mat')
            end
            if RMSEThRef ==3.5 
                    load('TheoreticalSpectra/A0_changingSNR.mat')
            end
            if RMSEThRef ==4 
                    load('TheoreticalSpectra/A0_changingSimilarity_400.mat')
            end
            if RMSEThRef ==4.5 
                    load('TheoreticalSpectra/A0_changingSimilarity_450.mat')
            end
            if RMSEThRef ==5 
                    load('TheoreticalSpectra/A0_changingSimilarity_500.mat')
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
    
    % Theta study
    if length(Vecttheta1)>1 | length(VectalphaA)>1
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
    
    % Data Preprocessing.
    [Y,Y_sps, H_sps] = dataPreprocessing(A0,Y);
    
    % Simplex projection of data.
    [UU, SS, VV] = svds(Y./sum(Y),size(A0,2));
    LowHaux = UU'*(Y./sum(Y));
    figure; plot(LowHaux(2,:),LowHaux(3,:),'b.'); hold on;%ALL pixels
    %LowPure = UU'*Y_purepixels; plot(LowPure(2,:),LowPure(3,:),'c.'); %Pure Pixels
    LowSparsest = UU'*(Y_sps./sum(Y_sps)); plot(LowSparsest(2,:),LowSparsest(3,:),'c.'); %Sparsest Pixels
    LowAnorm = UU'*Aref; hold on; plot(LowAnorm(2,:),LowAnorm(3,:),'go','MarkerSize',7); %Reference Spectra
    LowA0 = UU'*A0; hold on; plot(LowA0(2,:),LowA0(3,:),'k+','MarkerSize',7); %Initial Spectra
    xlabel('Subspace Vector 3'); ylabel('Subspace Vector 1'); title('Simplex Projection')        
for theta1=Vecttheta1
for theta2=Vecttheta2
for alphaA=VectalphaA
for alphaX=VectalphaX
for method = ChoosingMethod %To calculate Anmf.
        switch method             
            case 3 % Theoretical Spectra                
                A = A0;                                   
            case 5 % NMF-RI using theoretical Spectra                                  
                [A, H, iter] =  NMF_RI(Y_sps,A0,H_sps, theta1, theta2, alphaA, alphaX);                
            end
        
        % Sum-to-one normalization of every Blind spectral unmixing method.
        A = A*diag(1./(sum(A,1)+eps));

        % Calculate the resulted Hnmf, concentration matrix.
        Hnmf = max(100*eps,pinv(A)*Y); 
        
        % Calculate the resulted Href, concentration matrix.
        Href = max(100*eps,pinv(Aref)*Y);
        
        if length(VectnChannels)>1  
            Href = [GT,Hnmf(:,size(GT,2)+1:end)];
        end
        
        % Concentration matrix using the initialization spectra
        H0 = max(100*eps,pinv(A0)*Y);
                        
        % Creating concentration matrix, using the filter method.
        Hfilt = max(100*eps,pinv(AFilt)*Y);
        
        % Used to calculate 
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
        disp(['RelativeAccuracy using ',num2str(nc),' Channels: ',num2str(RACC(errorCounting)), ...
            'theta1:',num2str(theta1),'theta2:',num2str(theta2), 'alphaA', num2str(alphaA), 'alphaX', num2str(alphaX) ]);
        RACC_Names = [RACC_Names,{[currentfilename,'_Chnnl',num2str(nc),'_Spectra',Spectra]}];     
        switch method            
            case 1
                 RACC_NMFML = [RACC_NMFML, RACC(errorCounting)];                 
            case 2
                 RACC_NMFSB = [RACC_NMFSB, RACC(errorCounting)];   
            case 3
                 RACC_Theoretical = [RACC_Theoretical, RACC(errorCounting)];   
            case 4
                 RACC_NMFSB_Gauss = [RACC_NMFSB_Gauss, RACC(errorCounting)];   
            case 5
                 RACC_NMFRI = [RACC_NMFRI, RACC(errorCounting)];  
                 iterations_NMFRI = [iterations_NMFRI, iter];
            case 6 
                 RACC_NMFML_Exp = [RACC_NMFML_Exp, RACC(errorCounting)];   
            case 7 
                RACC_NMFRI_Exp = [RACC_NMFRI_Exp, RACC(errorCounting)];  
            case 8
                RACC_NMFSB_WithoutPixlSel = [RACC_NMFSB_WithoutPixlSel, RACC(errorCounting)];  
            case 9
                RACC_NMFML_WithoutPixlSel = [RACC_NMFML_WithoutPixlSel, RACC(errorCounting)];  
            case 10
                RACC_NMFRI_WithoutPixlSel = [RACC_NMFRI_WithoutPixlSel, RACC(errorCounting)];
            case 11
                RACC_GaussianInit = [RACC_GaussianInit, RACC(errorCounting)];
            case 12
                RACC_ExponentialInit = [RACC_ExponentialInit, RACC(errorCounting)];
            case 13 
                RACC_NMFTBCR = [RACC_NMFTBCR, RACC(errorCounting)]; 
            case 14
                RACC_NMFTBCR_Exp = [RACC_NMFTBCR_Exp, RACC(errorCounting)]; 
            case 15
                RACC_NMFTBCR_WithoutPixlSel = [RACC_NMFTBCR_WithoutPixlSel, RACC(errorCounting)]; 
            case 16
                RACC_NMFTBCR_GraphCut = [RACC_NMFTBCR_GraphCut, RACC(errorCounting)]; 
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
    case 'theta'
        RACC_NMFRI(1:36) = (RACC_NMFRI(1:36) + RACC_NMFRI(37:end))/2;
        figure;
        h = heatmap(Vecttheta1,Vecttheta2,reshape(RACC_NMFRI(1:36),[length(Vecttheta1),length(Vecttheta1)])); 
        xlabel('Theta 1 - IntraLayer stopping criteria')
        ylabel('Theta 2 - InterLayer stopping criteria')        
        colormap summer; caxis([0.5 0.95])

        iterations_NMFRI(1:36) = (iterations_NMFRI(1:36) + iterations_NMFRI(37:end))/2;
        figure; 
        h = heatmap(Vecttheta1,Vecttheta2,reshape(iterations_NMFRI(1:36),[length(Vecttheta1),length(Vecttheta1)])); 
        xlabel('Theta 1 - IntraLayer stopping criteria')
        ylabel('Theta 2 - InterLayer stopping criteria')        
        colormap(flipud(summer));
        caxis([0 100])
    case 'alpha'
%         RACC_NMFRI_reordered=zeros(length(VectalphaA)*length(VectalphaX),1);
%         for file=1:length(Vectnfiles) 
%             iter = 1;
%             for alphaA=1:length(VectalphaA)
%             for alphaX=1:length(VectalphaX)
%                 RACC_NMFRI_reordered(iter) = RACC_NMFRI_reordered(iter) + RACC_NMFRI(iter*file);
%                 iter = iter+1;
%             end
%             end
%         end
%         RACC_NMFRI_reordered=RACC_NMFRI_reordered/length(Vectnfiles);
%         
        RACC_NMFRI(1:36) = (RACC_NMFRI(1:36) + RACC_NMFRI(37:end))/2;
        figure;
        h = heatmap(VectalphaA,VectalphaX,reshape(RACC_NMFRI(1:36),[length(VectalphaA),length(VectalphaX)])); 
        xlabel('Alpha 2')
        ylabel('Alpha 1')        
        colormap(summer);
        caxis([0.5 0.95])
        
        iterations_NMFRI(1:36) = (iterations_NMFRI(1:36) + iterations_NMFRI(37:end))/2;
        figure; 
        h = heatmap(VectalphaA,VectalphaX,reshape(iterations_NMFRI(1:36),[length(VectalphaA),length(VectalphaX)])); 
        xlabel('Alpha 2')
        ylabel('Alpha 1')
        colormap(flipud(summer));
        caxis([0 100])
        zlabel('daf')
end  