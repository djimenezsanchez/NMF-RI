% This code unmixes flow cytometry data from AURORA cytometers.
% It follows NMF-RI as described in the paper "NMF-RI: Blind spectral unmixing of highly mixed multispectral flow and image cytometry"
% Flow cytometry data (.fcs) was converted to (.txt) using Infinicyt.
% We used emission and excitation data to create theoretical spectra.
% Written by Daniel Jimenez-Sanchez.
clear all;
close all;
clc;

%-------------------------------------------------------------%
% Folders to load data.
% Directory to load Single Color data, Y(Observation matrices).
SingleColorFolder = 'FlowCytometryData/LST_samples/Single color controls/Raw/';
% Directory to generate/load the theoretical Spectra
TheoreticalSpectraFolder = 'FlowCytometryData/LST_samples/TheoreticalSpectra/';
% Directory to load the mixed data, Y(Observation matrix).
SampleDataFolder = 'FlowCytometryData/LST_samples/Sample data/Raw/'; 
% Directory to save the unmixed data, H(Concentration matrix).
UnmixedDataFolder = 'FlowCytometryData/LST_samples/Sample data/Unmixed/'; 

sampleFiles = 2; % Number of samples to join for unmixing.
numChannels = 38; % Number of channels of AURORA Cytometer.
%-------------------------------------------------------------%

%% Creation of theoretical spectra using emission and excitation of the fluorochrome in certain wavelengths.

% Wavelength ranges for each spectral channel on the AURORA machine. 
load([TheoreticalSpectraFolder,'FiltWavelength.mat'])

% Load Emission and excitation.
EmissionFiles = dir([TheoreticalSpectraFolder,'*Em.txt']);
AbsorptionFiles = dir([TheoreticalSpectraFolder,'*Abs.txt']);
TheoreticalSpectra = zeros(size(EmissionFiles,1),size(FiltWavelength,2)); 

for i = 1:size(TheoreticalSpectra,1)
    EmissionFile = load([TheoreticalSpectraFolder,EmissionFiles(i).name]);
    AbsorptionFile = load([TheoreticalSpectraFolder,AbsorptionFiles(i).name]); 
    disp(EmissionFiles(i).name);
    for ii = 1:size(TheoreticalSpectra,2)
        
        % Calculate mean for emision/excitation with filters.
        TheoreticalSpectra(i,ii) = mean(EmissionFile((FiltWavelength(2,ii)>=EmissionFile(:,1)) & (EmissionFile(:,1)>=FiltWavelength(1,ii)) ,2));
        if isnan(TheoreticalSpectra(i,ii)); 
            TheoreticalSpectra(i,ii)=0; 
        end
        % Divide emmission to absorption.
        TheoreticalSpectra(i,ii) = TheoreticalSpectra(i,ii)*mean(AbsorptionFile((FiltWavelength(3,ii)>=AbsorptionFile(:,1)) & (AbsorptionFile(:,1)>=FiltWavelength(3,ii)) ,2));
        if isnan(TheoreticalSpectra(i,ii)); 
            TheoreticalSpectra(i,ii)=0; 
        end

    end
end

TheoreticalSpectra = TheoreticalSpectra'*diag(1./(sum(TheoreticalSpectra',1)+eps));
% figure; plot(TheoreticalSpectra);
save([TheoreticalSpectraFolder,'TheoreticalSpectra.mat'],'TheoreticalSpectra'),

%% Creation of Reference spectra using single-stained data.
% Load control/singlestained data.
SingleColorFiles = dir([SingleColorFolder,'*.txt']);
% Load Control Spectra created
load([SingleColorFolder,'ControlSpectraAll'],'ControlSpectraAll');
load([TheoreticalSpectraFolder,'TheoreticalSpectra.mat']);
TheoreticalSpectra(TheoreticalSpectra<0.00001) = 0.001;

% Open Autofluorescence file, to generate AF-Specra.
fileInfo = readtable([SingleColorFolder,SingleColorFiles(length(SingleColorFiles)).name],'Delimiter',';');
fileInfo_Array = table2array(fileInfo(:,1:numChannels));
ControlSpectra = nansum(fileInfo_Array)'; 
ControlSpectra = ControlSpectra'./sum(ControlSpectra',2);
ControlSpectraAll(:,length(SingleColorFiles)) =  ControlSpectra;

% Opening single-stained data.
for file=1:length(SingleColorFiles)-1

    % Reading .txt following the order -> first:a1...txt second:a2...txt , etc.
    fileInfo = readtable([SingleColorFolder,SingleColorFiles(file).name],'Delimiter',';');
    fileInfo_Array = table2array(fileInfo(:,1:numChannels));

    A0 = [TheoreticalSpectra(:,file),ControlSpectraAll(:,9)];
    H0=pinv(A0)*fileInfo_Array';
    % Eliminate AF from Single-color data
    [ControlSpectra, HnmfRI] = NMF_RI(fileInfo_Array',A0,H0);
    
    ControlSpectra = ControlSpectra'./sum(ControlSpectra',2);
    ControlSpectraAll(:,file) =  ControlSpectra(1,:);

    save([SingleColorFolder,SingleColorFiles(file).name(1:end-4),'_Converted'],'fileInfo_Array','ControlSpectra');
    disp([SingleColorFiles(file).name])
end

% figure; plot(ControlSpectraAll); title('All Spectra'); 
% xlabel('Spectral Channels'); ylabel('Normalized Intensity');
save([SingleColorFolder,'ControlSpectraAll'],'ControlSpectraAll');

%% Unmixing using Theoretical Spectra - NMF-RI
% Load Observation Matrix
ObservationFiles = dir([SampleDataFolder,'*data.txt']);
fileInfo_Array = [];
% Read 2 samples and join them to create observation Matrix, Y.
for i=1:sampleFiles
    fileInfo = readtable([SampleDataFolder,ObservationFiles(i).name],'Delimiter',';');
    fileInfo_Array = [fileInfo_Array;table2array(fileInfo(:,1:numChannels))];   
end
% Preprocessing step - Observation Matrix.
Y = [fileInfo_Array]';
Y(:,any(Y==4194304))=100*eps; % eliminate events that are saturated.
Y(Y<=0) = 100*eps; % negative values to zero.

% Load Control Spectra => Areference
load([SingleColorFolder,'ControlSpectraAll'],'ControlSpectraAll'); 
% Sum-to-one normalization per spectrum.
ControlSpectraAll = ControlSpectraAll*diag(1./(sum(ControlSpectraAll,1)+eps));

% Calculate Blood AF, it may be different from theoretical AF, because it
% was calculated using spheres. Blood AF and Sphere AF are different.
BloodAFFile = dir([SampleDataFolder,'Blood_AF.txt']);
fileInfo = readtable([SampleDataFolder,BloodAFFile(1).name],'Delimiter',';');
fileInfo_Array = table2array(fileInfo(:,1:numChannels));
AutoFluorescence = mean(fileInfo_Array);
AutoFluorescence = AutoFluorescence./sum(AutoFluorescence);
ControlSpectraAll(:,9) = AutoFluorescence'; 

% Load theoretical Spectra
load([TheoreticalSpectraFolder,'TheoreticalSpectra.mat']);
% Join theoretical spectra(spheres), with Blood AF. 
A0 = [TheoreticalSpectra, AutoFluorescence'];
% CAUTION: if a channel of A0 is equal to zero, NMF-RI would
% let it zero at the end of the optimization.
A0(A0<=0) = 0.0001;
A0 = A0*diag(1./(sum(A0,1)+eps));

% Data Preprocessing.
[Y_sps, H_sps] = dataPreprocessing(A0,Y);

% Data Preprocessing.
A0_straight = A0;
A0_straight(1:30,end) = 1/30; 
[Y_sps, H_sps] = dataPreprocessing(A0_straight,Y);
tic; [AnmfRI, ~] = NMF_RI_twohier(Y_sps,A0,H_sps); toc;

% Unmix data using NMF-RI
tic; [AnmfRI, ~] = NMF_RI(Y_sps,A0,H_sps); toc;
tic; [AnmfRI, ~] = NMF_RI_twohier(Y_sps,A0_straight,H_sps); toc;

displayResult(AnmfRI,ControlSpectraAll,A0_straight);
% displayResultBW(AnmfRI,ControlSpectraAll,A0);


%% Creating H to read the result using Infinicyt
AnmfRI=AnmfRI./sum(AnmfRI); % Normalize the spectra.
ControlSpectraAll=ControlSpectraAll./sum(ControlSpectraAll); % Normalize the spectra.

for i=1:sampleFiles
    % Which is the name of the file?
    NameFile = ObservationFiles(i).name;
    % Load Observation Matrix
    ObservationFiles = dir([SampleDataFolder,'*data.txt']);
    % Load Original 
    fileInfo = readtable([SampleDataFolder,ObservationFiles(i).name],'Delimiter',';');
    Y = table2array(fileInfo(:,1:38))';   
    % Events with intensity equal to 4194304 are saturated, wwe eliminate
    % them.
    Y(:,any(Y==4194304))=100*eps;
    % Set to zero negative values.
    Y(Y<=0) = 100*eps;
        
    % NMFRI Spectra
    H = pinv(AnmfRI)*Y; 
    % Apply percentile 99 to display
    H=(H./prctile(H,99,2)).*mean(Y(:));
    H_table = array2table(H','VariableNames',{'Pacific_BlueT','OC515','FITC','PE','PerCP_Cy5_5','PE_Cy7', ...
                                              'APC','APCC750','AF'});
    writetable(H_table,[UnmixedDataFolder,NameFile(1:end-4),'_NMFRIUnmixed','.txt'],'Delimiter',';');

    % Control Spectra
    H = pinv(ControlSpectraAll)*Y; 
    % Apply percentile 99 to display
    H=(H./prctile(H,99,2)).*mean(Y(:));
    H_table = array2table(H','VariableNames',{'Pacific_BlueT','OC515','FITC','PE','PerCP_Cy5_5','PE_Cy7', ...
                                              'APC','APCC750','AF'});
    writetable(H_table,[UnmixedDataFolder,NameFile(1:end-4),'_ControlUnmixed','.txt'],'Delimiter',';');
end


