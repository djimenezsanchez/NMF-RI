function [] = displayResult(AnmfRI,ControlSpectraAll,A0)

% First Unmixing - Display Curves while Unmixng - channels 1 to 16
figure; title('Spectra of UV Laser (405 nm)');
colors = distinguishable_colors(size(A0,2));
colors = [colors(1:3,:);colors(5:9,:)];
for i=1:size(AnmfRI,2)-1
        hold on; plot(ControlSpectraAll(1:16,i),'LineStyle','-','LineWidth',2,'Color',colors(i,:));
end
plot(ControlSpectraAll(1:16,i+1),'LineStyle','-','LineWidth',2,'Color',[0,0,0]);
ax = gca;
for i=1:size(A0,2)-1
    hold on; plot(A0(1:16,i),'LineStyle','--','LineWidth',1,'Color',colors(i,:));
end
plot(A0(1:16,i+1),'LineStyle','--','LineWidth',1,'Color',[0,0,0]);    
for i=1:size(A0,2)-1
    hold on; plot(AnmfRI(1:16,i),'LineStyle','-','Marker','o','LineWidth',1,'Color',colors(i,:));
end
plot(AnmfRI(1:16,i+1),'LineStyle','-','Marker','o','LineWidth',1,'Color',[0,0,0]);
ylabel('Intensity (Normalized)'); xlabel('Channels');                    
set(ax, 'XTick', 1:2:16)
set(ax, 'XTickLabel', [428,458,508,549,594,664,720,780])
ylabel('Intensity (Normalized)'); xlabel('Center wavelength of spectral channels (nm)');        
axis tight
box on;
% channels 17:30 - second laser
figure; title('Spectra of Blue Laser (488 nm)');
colors = distinguishable_colors(size(A0,2));
colors = [colors(1:3,:);colors(5:9,:)];
for i=1:size(A0,2)-1
    hold on; plot(ControlSpectraAll(17:30,i),'LineStyle','-','LineWidth',2,'Color',colors(i,:));
end
plot(ControlSpectraAll(17:30,i+1),'LineStyle','-','LineWidth',2,'Color',[0,0,0]);
 ax = gca;
for i=1:size(A0,2)-1
     hold on; plot(A0(17:30,i),'LineStyle','--','LineWidth',1,'Color',colors(i,:));
end
plot(A0(17:30,i+1),'LineStyle','-','LineWidth',1,'Color',[0,0,0]);    
for i=1:size(A0,2)-1
     hold on; plot(AnmfRI(17:30,i),'LineStyle','-','Marker','o','LineWidth',1,'Color',colors(i,:));
end
plot(AnmfRI(17:30,i+1),'LineStyle','-','Marker','o','LineWidth',1,'Color',[0,0,0]);
set(ax, 'XTick', 1:14)
set(ax, 'XTickLabel', [508,528,549,571,594,618,660,678,697,717,738,760,783,812])
ylabel('Intensity (Normalized)'); xlabel('Center wavelength of spectral channels (nm)');        
lgd.NumColumns = 3;
axis tight
box on
% channels 3138 - second laser
figure; title('Spectra of RED Laser (640 nm)');
colors = distinguishable_colors(size(A0,2));
colors = [colors(1:3,:);colors(5:9,:)];
for i=1:size(A0,2)-1
    hold on; plot(ControlSpectraAll(31:38,i),'LineStyle','-','LineWidth',2,'Color',colors(i,:));
end
plot(ControlSpectraAll(31:38,i+1),'LineStyle','-','LineWidth',2,'Color',[0,0,0]);
 ax = gca;
for i=1:size(A0,2)-1
     hold on; plot(A0(31:38,i),'LineStyle','--','LineWidth',1,'Color',colors(i,:));
end
plot(A0(31:38,i+1),'LineStyle','--','LineWidth',1,'Color',[0,0,0]);    
for i=1:size(A0,2)-1
     hold on; plot(AnmfRI(31:38,i),'LineStyle','-','Marker','o','LineWidth',1,'Color',colors(i,:));
end
plot(AnmfRI(31:38,i+1),'LineStyle','-','Marker','o','LineWidth',1,'Color',[0,0,0]);
set(ax, 'XTick', 1:8)
set(ax, 'XTickLabel', [660,678,697,717,738,760,783,812])
ylabel('Intensity (Normalized)'); xlabel('Center wavelength of spectral channels (nm)');        
    lgd = legend('Pacific BlueT Reference','OC515 Reference','FITC Reference','PE Reference','PerCP Cy5.5 Reference','PE Cy7 Reference', ...
                 'APC Reference','APC C750 Reference','AF Reference', ...
                 'Pacific BlueT Initialization','OC515 Initialization','FITC Initialization','PE Initialization','PerCP Cy5.5 Initialization','PE Cy7 Initialization', ...
                 'APC Initialization','APC C750 Initialization','AF Initialization', ...
                 'Pacific BlueT NMF-RI','OC515 NMF-RI','FITC NMF-RI','PE NMF-RI','PerCP Cy5.5 NMF-RI','PE Cy7 NMF-RI', ...
                 'APC NMF-RI','APC C750 NMF-RI','AF NMF-RI');    
    lgd.NumColumns = 3;
axis tight
box on

end