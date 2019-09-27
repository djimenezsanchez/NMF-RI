function [] = displayResultBW(AnmfRI,ControlSpectraAll,A0)

%% First Unmixing - Display Curves while Unmixng - channels 1 to 16
figure; 
hold on;
plot(ControlSpectraAll(1:16,1),'LineStyle','-','LineWidth',4,'Color',[0,0,0]);
plot(ControlSpectraAll(1:16,4),'LineStyle','-','LineWidth',4,'Color',[.4,.4,.4]);
plot(ControlSpectraAll(1:16,7),'LineStyle','-','LineWidth',4,'Color',[.8,.8,.8]);
ax = gca;
plot(A0(1:16,1),'LineStyle','--','LineWidth',3,'Color',[0,0,0]);
plot(A0(1:16,4),'LineStyle','--','LineWidth',3,'Color',[.4,.4,.4]);
plot(A0(1:16,7),'LineStyle','--','LineWidth',3,'Color',[.8,.8,.8]);

plot(AnmfRI(1:16,1),'LineStyle','-','Marker','o','LineWidth',3,'Color',[0,0,0]);
plot(AnmfRI(1:16,4),'LineStyle','-','Marker','o','LineWidth',3,'Color',[.4,.4,.4]);
plot(AnmfRI(1:16,7),'LineStyle','-','Marker','o','LineWidth',3,'Color',[.8,.8,.8]);

ylabel('Intensity (Normalized)'); xlabel('Channels');                    
set(ax, 'XTick', 1:2:16)
set(ax, 'XTickLabel', [428,458,508,549,594,664,720,780])
ylabel('Intensity (Normalized)'); xlabel('Center wavelength of spectral channels (nm)');  
lgd = legend('Pacific BlueT Reference','PE Reference','APC Reference', ...
                 'Pacific BlueT Initialization','PE Initialization','APC Initialization',...                 
                 'Pacific BlueT NMF-RI','PE NMF-RI','APC NMF-RI');   
lgd.NumColumns = 3;
axis tight
box on;

% First Unmixing - Display Curves while Unmixng - channels 1 to 16
figure; 
hold on;
plot(ControlSpectraAll(1:16,2),'LineStyle','-','LineWidth',4,'Color',[0,0,0]);
plot(ControlSpectraAll(1:16,5),'LineStyle','-','LineWidth',4,'Color',[.4,.4,.4]);
plot(ControlSpectraAll(1:16,8),'LineStyle','-','LineWidth',4,'Color',[.8,.8,.8]);
ax = gca;
plot(A0(1:16,2),'LineStyle','--','LineWidth',3,'Color',[0,0,0]);
plot(A0(1:16,5),'LineStyle','--','LineWidth',3,'Color',[.4,.4,.4]);
plot(A0(1:16,8),'LineStyle','--','LineWidth',3,'Color',[.8,.8,.8]);

plot(AnmfRI(1:16,2),'LineStyle','-','Marker','o','LineWidth',3,'Color',[0,0,0]);
plot(AnmfRI(1:16,5),'LineStyle','-','Marker','o','LineWidth',3,'Color',[.4,.4,.4]);
plot(AnmfRI(1:16,8),'LineStyle','-','Marker','o','LineWidth',3,'Color',[.8,.8,.8]);

ylabel('Intensity (Normalized)'); xlabel('Channels');                    
set(ax, 'XTick', 1:2:16)
set(ax, 'XTickLabel', [428,458,508,549,594,664,720,780])
ylabel('Intensity (Normalized)'); xlabel('Center wavelength of spectral channels (nm)');        
lgd = legend('OC515 Reference','PerCP Cy5.5 Reference','APC C750 Reference', ...
                 'OC515 Initialization','PerCP Cy5.5 Initialization','APC C750 Initialization',...                 
                 'OC515 NMF-RI','PerCP Cy5.5 NMF-RI','APC C750 NMF-RI');  
lgd.NumColumns = 3;
axis tight
box on;

% First Unmixing - Display Curves while Unmixng - channels 1 to 16
figure; 
hold on;
plot(ControlSpectraAll(1:16,3),'LineStyle','-','LineWidth',4,'Color',[0,0,0]);
plot(ControlSpectraAll(1:16,6),'LineStyle','-','LineWidth',4,'Color',[.4,.4,.4]);
plot(ControlSpectraAll(1:16,9),'LineStyle','-','LineWidth',4,'Color',[.8,.8,.8]);
ax = gca;
plot(A0(1:16,3),'LineStyle','--','LineWidth',3,'Color',[0,0,0]);
plot(A0(1:16,6),'LineStyle','--','LineWidth',3,'Color',[.4,.4,.4]);
plot(A0(1:16,9),'LineStyle','--','LineWidth',3,'Color',[.8,.8,.8]);
plot(AnmfRI(1:16,3),'LineStyle','-','Marker','o','LineWidth',3,'Color',[0,0,0]);
plot(AnmfRI(1:16,6),'LineStyle','-','Marker','o','LineWidth',3,'Color',[.4,.4,.4]);
plot(AnmfRI(1:16,9),'LineStyle','-','Marker','o','LineWidth',3,'Color',[.8,.8,.8]);

ylabel('Intensity (Normalized)'); xlabel('Channels');                    
set(ax, 'XTick', 1:2:16)
set(ax, 'XTickLabel', [428,458,508,549,594,664,720,780])
ylabel('Intensity (Normalized)'); xlabel('Center wavelength of spectral channels (nm)');        
lgd = legend('FITC Reference','PE Cy7 Reference','AF Reference', ...
                 'FITC Initialization','PE Cy7 Initialization','AF Initialization',...                 
                 'FITC NMF-RI','PE Cy7 NMF-RI','AF NMF-RI');   
lgd.NumColumns = 3;
axis tight
box on;

%% Second: channels 17:30 - second laser
figure; 
hold on;
plot(ControlSpectraAll(17:30,1),'LineStyle','-','LineWidth',4,'Color',[0,0,0]);
plot(ControlSpectraAll(17:30,4),'LineStyle','-','LineWidth',4,'Color',[.4,.4,.4]);
plot(ControlSpectraAll(17:30,7),'LineStyle','-','LineWidth',4,'Color',[.8,.8,.8]);
ax = gca;
plot(A0(17:30,1),'LineStyle','--','LineWidth',3,'Color',[0,0,0]);
plot(A0(17:30,4),'LineStyle','--','LineWidth',3,'Color',[.4,.4,.4]);
plot(A0(17:30,7),'LineStyle','--','LineWidth',3,'Color',[.8,.8,.8]);
plot(AnmfRI(17:30,1),'LineStyle','-','Marker','o','LineWidth',3,'Color',[0,0,0]);
plot(AnmfRI(17:30,4),'LineStyle','-','Marker','o','LineWidth',3,'Color',[.4,.4,.4]);
plot(AnmfRI(17:30,7),'LineStyle','-','Marker','o','LineWidth',3,'Color',[.8,.8,.8]);
ylabel('Intensity (Normalized)'); xlabel('Channels');                    
set(ax, 'XTick', 1:14)
set(ax, 'XTickLabel', [508,528,549,571,594,618,660,678,697,717,738,760,783,812])
ylabel('Intensity (Normalized)'); xlabel('Center wavelength of spectral channels (nm)');  
lgd = legend('Pacific BlueT Reference','PE Reference','APC Reference', ...
                 'Pacific BlueT Initialization','PE Initialization','APC Initialization',...                 
                 'Pacific BlueT NMF-RI','PE NMF-RI','APC NMF-RI'); 
lgd.NumColumns = 3;
axis tight
box on;

% Second: channels 17:30 - second laser
figure; 
hold on;
plot(ControlSpectraAll(17:30,2),'LineStyle','-','LineWidth',4,'Color',[0,0,0]);
plot(ControlSpectraAll(17:30,5),'LineStyle','-','LineWidth',4,'Color',[.4,.4,.4]);
plot(ControlSpectraAll(17:30,8),'LineStyle','-','LineWidth',4,'Color',[.8,.8,.8]);
ax = gca;
plot(A0(17:30,2),'LineStyle','--','LineWidth',3,'Color',[0,0,0]);
plot(A0(17:30,5),'LineStyle','--','LineWidth',3,'Color',[.4,.4,.4]);
plot(A0(17:30,8),'LineStyle','--','LineWidth',3,'Color',[.8,.8,.8]);
plot(AnmfRI(17:30,2),'LineStyle','-','Marker','o','LineWidth',3,'Color',[0,0,0]);
plot(AnmfRI(17:30,5),'LineStyle','-','Marker','o','LineWidth',3,'Color',[.4,.4,.4]);
plot(AnmfRI(17:30,8),'LineStyle','-','Marker','o','LineWidth',3,'Color',[.8,.8,.8]);
ylabel('Intensity (Normalized)'); xlabel('Channels');                    
set(ax, 'XTick', 1:14)
set(ax, 'XTickLabel', [508,528,549,571,594,618,660,678,697,717,738,760,783,812])
ylabel('Intensity (Normalized)'); xlabel('Center wavelength of spectral channels (nm)');        
lgd = legend('OC515 Reference','PerCP Cy5.5 Reference','APC C750 Reference', ...
                 'OC515 Initialization','PerCP Cy5.5 Initialization','APC C750 Initialization',...                 
                 'OC515 NMF-RI','PerCP Cy5.5 NMF-RI','APC C750 NMF-RI');
lgd.NumColumns = 3;
axis tight
box on;

% Second: channels 17:30 - second laser
figure; 
hold on;
plot(ControlSpectraAll(17:30,3),'LineStyle','-','LineWidth',4,'Color',[0,0,0]);
plot(ControlSpectraAll(17:30,6),'LineStyle','-','LineWidth',4,'Color',[.4,.4,.4]);
plot(ControlSpectraAll(17:30,9),'LineStyle','-','LineWidth',4,'Color',[.8,.8,.8]);
ax = gca;
plot(A0(17:30,3),'LineStyle','--','LineWidth',3,'Color',[0,0,0]);
plot(A0(17:30,6),'LineStyle','--','LineWidth',3,'Color',[.4,.4,.4]);
plot(A0(17:30,9),'LineStyle','--','LineWidth',3,'Color',[.8,.8,.8]);
plot(AnmfRI(17:30,3),'LineStyle','-','Marker','o','LineWidth',3,'Color',[0,0,0]);
plot(AnmfRI(17:30,6),'LineStyle','-','Marker','o','LineWidth',3,'Color',[.4,.4,.4]);
plot(AnmfRI(17:30,9),'LineStyle','-','Marker','o','LineWidth',3,'Color',[.8,.8,.8]);
ylabel('Intensity (Normalized)'); xlabel('Channels');                    
set(ax, 'XTick', 1:14)
set(ax, 'XTickLabel', [508,528,549,571,594,618,660,678,697,717,738,760,783,812])
ylabel('Intensity (Normalized)'); xlabel('Center wavelength of spectral channels (nm)');        
lgd = legend('FITC Reference','PE Cy7 Reference','AF Reference', ...
                 'FITC Initialization','PE Cy7 Initialization','AF Initialization',...                 
                 'FITC NMF-RI','PE Cy7 NMF-RI','AF NMF-RI');    
lgd.NumColumns = 3;
axis tight
box on;

%% channels 3138 - third laser
figure; 
hold on;
plot(ControlSpectraAll(31:38,1),'LineStyle','-','LineWidth',4,'Color',[0,0,0]);
plot(ControlSpectraAll(31:38,4),'LineStyle','-','LineWidth',4,'Color',[.4,.4,.4]);
plot(ControlSpectraAll(31:38,7),'LineStyle','-','LineWidth',4,'Color',[.8,.8,.8]);
ax = gca;
plot(A0(31:38,1),'LineStyle','--','LineWidth',3,'Color',[0,0,0]);
plot(A0(31:38,4),'LineStyle','--','LineWidth',3,'Color',[.4,.4,.4]);
plot(A0(31:38,7),'LineStyle','--','LineWidth',3,'Color',[.8,.8,.8]);
plot(AnmfRI(31:38,1),'LineStyle','-','Marker','o','LineWidth',3,'Color',[0,0,0]);
plot(AnmfRI(31:38,4),'LineStyle','-','Marker','o','LineWidth',3,'Color',[.4,.4,.4]);
plot(AnmfRI(31:38,7),'LineStyle','-','Marker','o','LineWidth',3,'Color',[.8,.8,.8]);
ylabel('Intensity (Normalized)'); xlabel('Channels');                    
set(ax, 'XTick', 1:8)
set(ax, 'XTickLabel', [660,678,697,717,738,760,783,812])
ylabel('Intensity (Normalized)'); xlabel('Center wavelength of spectral channels (nm)');        
lgd = legend('Pacific BlueT Reference','PE Reference','APC Reference', ...
                 'Pacific BlueT Initialization','PE Initialization','APC Initialization',...                 
                 'Pacific BlueT NMF-RI','PE NMF-RI','APC NMF-RI');    
lgd.NumColumns = 3;
axis tight
box on;

% Second: channels 17:30 - second laser
figure; 
hold on;
plot(ControlSpectraAll(31:38,2),'LineStyle','-','LineWidth',4,'Color',[0,0,0]);
plot(ControlSpectraAll(31:38,5),'LineStyle','-','LineWidth',4,'Color',[.4,.4,.4]);
plot(ControlSpectraAll(31:38,8),'LineStyle','-','LineWidth',4,'Color',[.8,.8,.8]);
ax = gca;
plot(A0(31:38,2),'LineStyle','--','LineWidth',3,'Color',[0,0,0]);
plot(A0(31:38,5),'LineStyle','--','LineWidth',3,'Color',[.4,.4,.4]);
plot(A0(31:38,8),'LineStyle','--','LineWidth',3,'Color',[.8,.8,.8]);
plot(AnmfRI(31:38,2),'LineStyle','-','Marker','o','LineWidth',3,'Color',[0,0,0]);
plot(AnmfRI(31:38,5),'LineStyle','-','Marker','o','LineWidth',3,'Color',[.4,.4,.4]);
plot(AnmfRI(31:38,8),'LineStyle','-','Marker','o','LineWidth',3,'Color',[.8,.8,.8]);
ylabel('Intensity (Normalized)'); xlabel('Channels');                    
set(ax, 'XTick', 1:8)
set(ax, 'XTickLabel', [660,678,697,717,738,760,783,812])
ylabel('Intensity (Normalized)'); xlabel('Center wavelength of spectral channels (nm)');        
lgd = legend('OC515 Reference','PerCP Cy5.5 Reference','APC C750 Reference', ...
                 'OC515 Initialization','PerCP Cy5.5 Initialization','APC C750 Initialization',...                 
                 'OC515 NMF-RI','PerCP Cy5.5 NMF-RI','APC C750 NMF-RI');
lgd.NumColumns = 3;
axis tight
box on;

% Second: channels 17:30 - second laser
figure; 
hold on;
plot(ControlSpectraAll(31:38,3),'LineStyle','-','LineWidth',4,'Color',[0,0,0]);
plot(ControlSpectraAll(31:38,6),'LineStyle','-','LineWidth',4,'Color',[.4,.4,.4]);
plot(ControlSpectraAll(31:38,9),'LineStyle','-','LineWidth',4,'Color',[.8,.8,.8]);
ax = gca;
plot(A0(31:38,3),'LineStyle','--','LineWidth',3,'Color',[0,0,0]);
plot(A0(31:38,6),'LineStyle','--','LineWidth',3,'Color',[.4,.4,.4]);
plot(A0(31:38,9),'LineStyle','--','LineWidth',3,'Color',[.8,.8,.8]);
plot(AnmfRI(31:38,3),'LineStyle','-','Marker','o','LineWidth',3,'Color',[0,0,0]);
plot(AnmfRI(31:38,6),'LineStyle','-','Marker','o','LineWidth',3,'Color',[.4,.4,.4]);
plot(AnmfRI(31:38,9),'LineStyle','-','Marker','o','LineWidth',3,'Color',[.8,.8,.8]);
ylabel('Intensity (Normalized)'); xlabel('Channels');                    
set(ax, 'XTick', 1:8)
set(ax, 'XTickLabel', [660,678,697,717,738,760,783,812])
ylabel('Intensity (Normalized)'); xlabel('Center wavelength of spectral channels (nm)');        
lgd = legend('FITC Reference','PE Cy7 Reference','AF Reference', ...
                 'FITC Initialization','PE Cy7 Initialization','AF Initialization',...                 
                 'FITC NMF-RI','PE Cy7 NMF-RI','AF NMF-RI');    
lgd.NumColumns = 3;
axis tight
box on;

end