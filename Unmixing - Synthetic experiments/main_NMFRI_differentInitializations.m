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
        ChoosingMethod = [1,2,3,4,5,6,7,8,9,10]; % Number of methods to compare        
        %%%%%%%%%%%%%%%%%%%%%%%-Parameters-%%%%%%%%%%%%%%%%%%%%%%%
    case 'RMSEThRef'
        %%%%%%%%%%%%%%%%%%%%%%%-Parameters-%%%%%%%%%%%%%%%%%%%%%%%
        % These parameters studies different RMSE between Theoretical and Reference Spectra.
        nf = 3;                     % number of fluorochromes        
        VectSNR = 35;               % SNR levels
        VectRMSEThRef = [4,4.25, 4.5];    % Similarity between Reference and theoretical spectra
        Vectnfiles = [4];             % Files to open with different levels of colocalization
        VectnChannels=36;           % Number of channels              
        ChoosingMethod = [1,2,3]; % Number of methods to compare        
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

RACC_Theoretical=[]; Time_Theoretical=[];
RACC_NMFRI=[]; Time_NMFRI=[];
RACC_NMFRI_Exp=[]; Time_NMFRI_Exp=[];
RACC_NMFRI_Gauss=[]; Time_NMFRI_Gauss=[];
RACC_NMFRI_GraphCut=[]; Time_NMFRI_GraphCut=[];
RACC_NMFRI_VCA=[]; Time_NMFRI_VCA=[];
RACC_Names=[]; elapsedTime=[];
RACC_Exp=[]; Time_Exp=[];
RACC_Gauss=[]; Time_Gauss=[];
RACC_GraphCut=[]; Time_GraphCut=[];
RACC_random=[]; Time_random=[];
RACC_Initrandom=[]; Time_Initrandom=[];
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
            case 1 % Theoretical Spectra                
                A = A0;                       
            case 2 % NMF-RI using theoretical Spectra                                  
                [A, H, elapsedTime] =  NMF_RI(Y_sps,A0,H_sps, theta1, theta2, alphaA, alphaX);            
            case 3 % NMF-RI using Exponential Spectra
                % In this case we recalculate the H
                H0exp = max(100*eps,pinv(exponential_matrix(nf,nc))*Y_sps); 
                [A, H,elapsedTime] = NMF_RI(Y_sps,exponential_matrix(nf,nc),H0exp, theta1, theta2, alphaA, alphaX);
            case 4 % NMF-RI using Gaussian Spectra
                H0gauss = max(100*eps,pinv(gauss_matrix(nf,nc))*Y_sps); 
                [A, H,elapsedTime] =  NMF_RI(Y_sps,gauss_matrix(nf,nc),H0gauss, theta1, theta2, alphaA, alphaX);                   
            case 5 % NMF-RI with Graph-Cut pixel selection 
                addpath('TBCR-NMF')
                ACutInit = graphCutInit(Y_sps,size(A0,2));                
                Hinit = max(100*eps,pinv(ACutInit)*Y_sps); 
                [A, H, elapsedTime] =  NMF_RI(Y_sps,ACutInit,Hinit, theta1, theta2, alphaA, alphaX);   
            case 6 % Exponential spectra
                A = exponential_matrix(nf,nc);
            case 7 % Gaussian spectra
                A = gauss_matrix(nf,nc);
            case 8 % nCUT spectra
                A = ACutInit;
            case 9 % Random spectra
                Aend = Aref*0;
                for iii = 1:50
                    A0rand = rand(size(Aref));
                    A0rand = A0rand./sum(A0rand,1);
                    Hinit = max(100*eps,pinv(A0rand)*Y_sps); 
                    [A, H, elapsedTime] =  NMF_RI(Y_sps,A0rand,Hinit, 0.001, 0.001, alphaA, 0.1);   
                    [D] = pdist2(Aref',A','euclidean');
                    [~, D] = min(D,[],2);
                    A=A(:,D);
                    Aend = A+Aend;
                end
                A = Aend;
            case 10 % Random initialization
                Aend = Aref*0;
                for iii = 1:50
                    A0rand = rand(size(Aref));
                    A0rand = A0rand./sum(A0rand,1);
                    [D] = pdist2(Aref',A0rand','euclidean');
                    [~, D] = min(D,[],2);
                    A0rand=A0rand(:,D);
                    Aend = A0rand+Aend;
                end
                A = Aend;
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
        disp(['RelativeAccuracy using ',num2str(nc),' Channels: ',num2str(RACC(errorCounting) )]);
        RACC_Names = [RACC_Names,{[currentfilename,'_Chnnl',num2str(nc),'_Spectra',Spectra]}];
        % Values that are negative put to zero.
%         RACC(RACC<0)=0;
        switch method            
            case 1
                RACC_Theoretical = [RACC_Theoretical, RACC(errorCounting)];   
                 Time_Theoretical = [Time_Theoretical, elapsedTime];                      
            case 2
                RACC_NMFRI = [RACC_NMFRI, RACC(errorCounting)];                 
                 Time_NMFRI = [Time_NMFRI, elapsedTime];              
            case 3
                 RACC_NMFRI_Exp = [RACC_NMFRI_Exp, RACC(errorCounting)];   
                 Time_NMFRI_Exp = [Time_NMFRI_Exp, elapsedTime];   
                 
            case 4
                RACC_NMFRI_Gauss = [RACC_NMFRI_Gauss, RACC(errorCounting)];   
                 Time_NMFRI_Gauss = [Time_NMFRI_Gauss, elapsedTime];   
                    
            case 5
                RACC_NMFRI_GraphCut = [RACC_NMFRI_GraphCut, RACC(errorCounting)];   
                 Time_NMFRI_GraphCut = [Time_NMFRI_GraphCut, elapsedTime];   
            case 6
                RACC_Exp = [RACC_Exp, RACC(errorCounting)];   
                 Time_Exp = [Time_Exp, elapsedTime];   
            case 7
                RACC_Gauss = [RACC_Gauss, RACC(errorCounting)];   
                Time_Gauss = [Time_Gauss, elapsedTime];   
            case 8
                RACC_GraphCut = [RACC_GraphCut, RACC(errorCounting)];   
                Time_GraphCut = [Time_GraphCut, elapsedTime];
            case 9 
                RACC_random = [RACC_random, RACC(errorCounting)];   
                Time_random = [Time_random, elapsedTime];
            case 10
                RACC_Initrandom = [RACC_Initrandom, RACC(errorCounting)];   
                Time_Initrandom = [Time_Initrandom, elapsedTime];
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
exp = exponential_matrix(nf,nc);
exp = exp./sum(exp);
gaus = gauss_matrix(nf,nc);
gaus = gaus./sum(gaus);
plot(vectorA,Aref(:,2),'Color','b','LineStyle','-','LineWidth',1);
plot(vectorA,Aref(:,3),'Color','r','LineStyle','-','LineWidth',1);
plot(vectorA,exp(:,1),'Color','g','Marker','o','LineStyle','--','LineWidth',1);
plot(vectorA,exp(:,2),'Color','b','Marker','o','LineStyle','--','LineWidth',1);
plot(vectorA,exp(:,3),'Color','r','Marker','o','LineStyle','--','LineWidth',1);
plot(vectorA,gaus(:,1),'Color','g','Marker','*','LineStyle','--','LineWidth',1);
plot(vectorA,gaus(:,2),'Color','b','Marker','*','LineStyle','--','LineWidth',1);
plot(vectorA,gaus(:,3),'Color','r','Marker','*','LineStyle','--','LineWidth',1);
plot(vectorA,A0(:,1),'Color','g','LineStyle','--','LineWidth',1);
plot(vectorA,A0(:,2),'Color','b','LineStyle','--','LineWidth',1);
plot(vectorA,A0(:,3),'Color','r','LineStyle','--','LineWidth',1);
legend('Alexa488 Reference','R-PE Reference','PE-cy5 Reference','Alexa488 Exponential Initialization',...
    'R-PE Exponential Initialization','PE-cy5 Exponential Initialization','Alexa488 Gaussian Initialization',...
    'R-PE Gaussian Initialization','PE-cy5 Gaussian Initialization', ...
    'Alexa488 Theoretical Initialization','R-PE Theoretical Initialization','PE-cy5 Theoretical Initialization');
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
        plot(ColLevels,RACC_NMFRI_Gauss,'LineWidth',2,'LineStyle','-'); 
        plot(ColLevels,RACC_NMFRI_Exp,'LineWidth',2,'LineStyle','-'); 
        plot(ColLevels,RACC_NMFRI,'LineWidth',2,'LineStyle','-');
        plot(ColLevels,RACC_NMFRI_GraphCut,'LineWidth',2,'LineStyle','-');  
        plot(ColLevels,RACC_random,'LineWidth',2,'LineStyle','-','Color',[0, 0, 0]);                                   
        plot(ColLevels,RACC_Gauss,'LineWidth',2,'LineStyle','--','Color',[0.4940, 0.1840, 0.5560]);
        plot(ColLevels,RACC_Exp,'LineWidth',2,'LineStyle','--','Color',[0.4660, 0.6740, 0.1880]);
        plot(ColLevels,RACC_GraphCut,'LineWidth',2,'LineStyle','--','Color',[0.6350, 0.0780, 0.1840]);
        plot(ColLevels,RACC_Initrandom,'LineWidth',2,'LineStyle','--','Color',[0, 0, 0]);     
        legend('Reference Spectra','Filter method','Theoretical Spectra','NMF-RI with Gauss Initialization','NMF-RI with Exponential Initialization','NMF-RI withTheoretical Initialization',...
            'NMF-RI with nCUT Initialization', 'NMF-RI with random initialization', 'Gauss initialization spectra', 'Exponential initialization spectra', 'nCUT initialization spectra', 'Random initialization spectra');
        xlabel('Degrees of spatial colocalization (%)');
        ylabel('RACC');
        xticks([0,25,50,75,100]); xticklabels({'0','25','50','75','100'});
 
end  