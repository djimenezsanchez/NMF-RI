% This code unmixes Image Cytometry data from Vectra Polaris system.
% It follows NMF-RI as described in the paper "NMF-RI: Blind spectral unmixing of highly mixed multispectral flow and image cytometry"
% Observation matrices were obtained from Vectra Polaris as a tiff-stack.
% We used emission and excitation data to create theoretical spectra.
% Written by Daniel Jimenez-Sanchez.

clear all;
close all;
clc;

%-------------------------------------------------------------%
% Folders to load data.
% Directory to load Single Color data, Y(Observation matrices).
SingleColorFolder = 'HumanOvary/Single color controls/';
% Directory to generate/load the theoretical Spectra
TheoreticalSpectraFolder = 'HumanOvary/TheoreticalSpectra/';
% Directory to load the mixed data from .tiff, Y(Observation matrix).
SampleDataFolderOriginal = 'HumanOvary/Sample data/raw/OriginalImages/';
% Directory to load the mixed data from .mat, Y(Observation matrix).
SampleDataFolder = 'HumanOvary/Sample data/raw/'; 
% Directory to save the unmixed data, H(Concentration matrix).
UnmixedDataFolder = 'HumanOvary/Sample data/unmixed/'; 

sampleFiles = 5; % Number of images to join for unmixing.
%-------------------------------------------------------------%

%% Creating Observation Matrix.
Testfiles  = dir([SampleDataFolderOriginal,'*.tif']);

chnl=1;
for file=1:length(Testfiles)
for i=1:9 % First Filter
IM_1(:,:,chnl) = imread([SampleDataFolderOriginal,Testfiles(file).name],i); chnl=chnl+1;
end
IM_1=IM_1;
[IM_1, IM_1_MASK] = toMask(IM_1); chnl=1; 
for i=10:18 % Second Filter
IM_2(:,:,chnl) = imread([SampleDataFolderOriginal,Testfiles(file).name],i); chnl=chnl+1;
end
IM_2=IM_2;
[IM_2, IM_2_MASK] = toMask(IM_2); chnl=1;
for i=19:25 % Third Filter
IM_3(:,:,chnl) = imread([SampleDataFolderOriginal,Testfiles(file).name],i); chnl=chnl+1;
end
IM_3=IM_3;
[IM_3, IM_3_MASK] = toMask(IM_3); chnl=1;
for i=26:32 % Fourth Filter
IM_4(:,:,chnl) = imread([SampleDataFolderOriginal,Testfiles(file).name],i); chnl=chnl+1;
end
IM_4=IM_4;
[IM_4, IM_4_MASK] = toMask(IM_4); chnl=1;
for i=33:35 % Fifth Filter
IM_5(:,:,chnl) = imread([SampleDataFolderOriginal,Testfiles(file).name],i); chnl=chnl+1;
end
IM_5=IM_5;
[IM_5, IM_5_MASK] = toMask(IM_5); chnl=1;

save([SampleDataFolder,Testfiles(file).name(1:end-4)],'IM_1','IM_1_MASK','IM_2','IM_2_MASK','IM_3','IM_3_MASK','IM_4','IM_4_MASK','IM_5','IM_5_MASK');
clear IM_1 IM_1_MASK IM_2 IM_2_MASK IM_3 IM_3_MASK IM_4 IM_4_MASK IM_5 IM_5_MASK;
disp(Testfiles(file).name)
end


%% Unmixing using Theoretical Spectra - NMF-RI
% Load Observation Matrix
ObservationFiles = dir([SampleDataFolder,'*.mat']);
Y = [];
% Read 5 images and join them to create an observation Matrix, Y.
for iImage=1:sampleFiles    
    ObsFile = load([SampleDataFolder,ObservationFiles(iImage).name]);
    Ynext = [ObsFile.IM_1,ObsFile.IM_2,ObsFile.IM_3,ObsFile.IM_4,ObsFile.IM_5]';
    Y = [Y, Ynext];
    clear Ynext
end
Y = Y-medfilt1(min(Y));
Y(Y<=0) = 100*eps; % negative values to zero.
    
% Load ControlSpectraAll => Areference. It was obtained from single-stained
% samples prcessed using Inform software. (PerkinElmer, CA, USA.)
load([SingleColorFolder,'ControlSpectra']); 
% Sum-to-one normalization per spectrum.
ControlSpectraAll = ControlSpectraAll*diag(1./(sum(ControlSpectraAll,1)+eps));

% Load theoretical Spectra
load([TheoreticalSpectraFolder,'TheoreticalSpectra.mat']);
TheoreticalSpectra=TheoreticalSpectra*diag(1./(sum(TheoreticalSpectra,1)+eps));
% Join theoretical spectra(spheres), with Control AF. 
A0 = [TheoreticalSpectra(:,1:7), ControlSpectraAll(:,8)];

% Data Preprocessing.
[Y_sps, H_sps] = dataPreprocessing(A0,Y);
   
% Unmix data using NMF-RI
tic; [AnmfRI, ~] = NMF_RI(Y_sps,A0,H_sps); toc;

% Display the results of each filter
displayResult(AnmfRI,ControlSpectraAll,A0);

%% Calculate and save result.

for iImage=1:sampleFiles    
    % Load Observation Matrix
    ObservationFiles = dir([SampleDataFolder,'*.mat']);
    ObsFile = load([SampleDataFolder,ObservationFiles(iImage).name]);
    Y = [ObsFile.IM_1,ObsFile.IM_2,ObsFile.IM_3,ObsFile.IM_4,ObsFile.IM_5]';
    
    % Linear unmixing. NMF-RI
    H_nmri = max(1E6*eps,pinv(AnmfRI'*AnmfRI)*AnmfRI'*Y);  
    % zscore normalization
    H_nmri = zscore(H_nmri')'-min(zscore(H_nmri')');
    
    % Linear unmixing. Controls
    H_control = max(1E6*eps,pinv(ControlSpectraAll'*ControlSpectraAll)*ControlSpectraAll'*Y);   
    % zscore normalization
    H_control = zscore(H_control')'-min(zscore(H_control')');

H_nmri = (H_nmri./max(H_nmri,[],2))*65535; H_control = (H_control./max(H_control,[],2))*65535;
IM__nmfri = uint16(double(reshape(H_nmri',[1404,1876,8]))); % Reshape to the Size of the original images
IM_controls = uint16(double(reshape(H_control',[1404,1876,8])));

% Write the images (H) to the UnmixedDataFolder directory
imwrite(IM__nmfri(:,:,1), [UnmixedDataFolder,'Ovary_NMF_Unmixed',num2str(iImage),'.TIFF'])
imwrite(IM__nmfri(:,:,2), [UnmixedDataFolder,'Ovary_NMF_Unmixed',num2str(iImage),'.TIFF'],'writemode', 'append')
imwrite(IM__nmfri(:,:,3), [UnmixedDataFolder,'Ovary_NMF_Unmixed',num2str(iImage),'.TIFF'], 'writemode', 'append')
imwrite(IM__nmfri(:,:,4), [UnmixedDataFolder,'Ovary_NMF_Unmixed',num2str(iImage),'.TIFF'], 'writemode', 'append')
imwrite(IM__nmfri(:,:,5), [UnmixedDataFolder,'Ovary_NMF_Unmixed',num2str(iImage),'.TIFF'], 'writemode', 'append')
imwrite(IM__nmfri(:,:,6), [UnmixedDataFolder,'Ovary_NMF_Unmixed',num2str(iImage),'.TIFF'], 'writemode', 'append')
imwrite(IM__nmfri(:,:,7), [UnmixedDataFolder,'Ovary_NMF_Unmixed',num2str(iImage),'.TIFF'], 'writemode', 'append')

% Write the images (H) to the UnmixedDataFolder directory
imwrite(IM_controls(:,:,1), [UnmixedDataFolder,'Ovary_Controls_Unmixed',num2str(iImage),'.TIFF'])
imwrite(IM_controls(:,:,2), [UnmixedDataFolder,'Ovary_Controls_Unmixed',num2str(iImage),'.TIFF'], 'writemode', 'append')
imwrite(IM_controls(:,:,3), [UnmixedDataFolder,'Ovary_Controls_Unmixed',num2str(iImage),'.TIFF'], 'writemode', 'append')
imwrite(IM_controls(:,:,4), [UnmixedDataFolder,'Ovary_Controls_Unmixed',num2str(iImage),'.TIFF'], 'writemode', 'append')
imwrite(IM_controls(:,:,5), [UnmixedDataFolder,'Ovary_Controls_Unmixed',num2str(iImage),'.TIFF'], 'writemode', 'append')
imwrite(IM_controls(:,:,6), [UnmixedDataFolder,'Ovary_Controls_Unmixed',num2str(iImage),'.TIFF'], 'writemode', 'append')
imwrite(IM_controls(:,:,7), [UnmixedDataFolder,'Ovary_Controls_Unmixed',num2str(iImage),'.TIFF'], 'writemode', 'append')    

end
