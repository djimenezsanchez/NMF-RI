function [] = displayResult(AnmfRI,ControlSpectraAll,A0)

%% figure 1
figure; title('Spectra on DAPI Filter');  ax=gca;
hold on; plot(ControlSpectraAll(1:9,1),'LineStyle','-','LineWidth',2,'Color',[0,0,1]);
hold on; plot(ControlSpectraAll(1:9,2),'LineStyle','-','LineWidth',2,'Color',[0,1,0]);
hold on; plot(ControlSpectraAll(1:9,3),'LineStyle','-','LineWidth',2,'Color',[.9,.9,0]);
hold on; plot(ControlSpectraAll(1:9,4),'LineStyle','-','LineWidth',2,'Color',[1,0,0]);
hold on; plot(ControlSpectraAll(1:9,5),'LineStyle','-','LineWidth',2,'Color',[1,0,1]);
hold on; plot(ControlSpectraAll(1:9,6),'LineStyle','-','LineWidth',2,'Color',[.5,.5,.5]);
hold on; plot(ControlSpectraAll(1:9,7),'LineStyle','-','LineWidth',2,'Color',[0,1,1]);
hold on; plot(ControlSpectraAll(1:9,8),'LineStyle','-','LineWidth',2,'Color',[0,0,0]);

hold on; plot(A0(1:9,1),'LineStyle','--','LineWidth',1,'Color',[0,0,1]);
hold on; plot(A0(1:9,2),'LineStyle','--','LineWidth',1,'Color',[0,1,0]);
hold on; plot(A0(1:9,3),'LineStyle','--','LineWidth',1,'Color',[.9,.9,0]);
hold on; plot(A0(1:9,4),'LineStyle','--','LineWidth',1,'Color',[1,0,0]);
hold on; plot(A0(1:9,5),'LineStyle','--','LineWidth',1,'Color',[1,0,1]);
hold on; plot(A0(1:9,6),'LineStyle','--','LineWidth',1,'Color',[.5,.5,.5]);
hold on; plot(A0(1:9,7),'LineStyle','--','LineWidth',1,'Color',[0,1,1]);
hold on; plot(A0(1:9,8),'LineStyle','--','LineWidth',1,'Color',[0,0,0]);

hold on; plot(AnmfRI(1:9,1),'LineStyle','-','Marker','*','LineWidth',1,'Color',[0,0,1]);
hold on; plot(AnmfRI(1:9,2),'LineStyle','-','Marker','*','LineWidth',1,'Color',[0,1,0]);
hold on; plot(AnmfRI(1:9,3),'LineStyle','-','Marker','*','LineWidth',1,'Color',[.9,.9,0]);
hold on; plot(AnmfRI(1:9,4),'LineStyle','-','Marker','*','LineWidth',1,'Color',[1,0,0]);
hold on; plot(AnmfRI(1:9,5),'LineStyle','-','Marker','*','LineWidth',1,'Color',[1,0,1]);
hold on; plot(AnmfRI(1:9,6),'LineStyle','-','Marker','*','LineWidth',1,'Color',[.5,.5,.5]);
hold on; plot(AnmfRI(1:9,7),'LineStyle','-','Marker','*','LineWidth',1,'Color',[0,1,1]);
hold on; plot(AnmfRI(1:9,8),'LineStyle','-','Marker','*','LineWidth',1,'Color',[0,0,0]);
set(ax,'XTick',1:9);
set(ax,'XTickLabel',[440,460,480,500,520,540,560,580,600])
xlabel('Center wavelength of spectral channels (nm)')
ylabel('Intensity (Normalized)')
axis tight;
box on;

h = legend('DAPI Reference', 'Opal 520 Reference', 'Opal 540 Reference', 'Opal 570 Reference', ...
    'Opal 620 Reference','Opal 650 Reference', 'Opal 690 Reference', 'AF Reference', ...
    'DAPI Initialization', 'Opal 520 Initialization', 'Opal 540 Initialization', 'Opal 570 Initialization', ...
    'Opal 620 Initialization','Opal 650 Initialization', 'Opal 690 Initialization', 'AF Initialization', ...
    'DAPI NMF-RI (Without AF reference)', 'Opal 520 NMF-RI (Without AF reference)', ...
    'Opal 540 NMF-RI (Without AF reference)', 'Opal 570 NMF-RI (Without AF reference)', ...
    'Opal 620 NMF-RI (Without AF reference)','Opal 650 NMF-RI (Without AF reference)',...
    'Opal 690 NMF-RI (Without AF reference)', 'AF NMF-RI (Without AF reference)' ...
    );
h.NumColumns = 3;


%% figure 2
figure; title('Spectra on FITC Filter');  ax=gca;
hold on; plot(ControlSpectraAll(10:18,1),'LineStyle','-','LineWidth',2,'Color',[0,0,1]);
hold on; plot(ControlSpectraAll(10:18,2),'LineStyle','-','LineWidth',2,'Color',[0,1,0]);
hold on; plot(ControlSpectraAll(10:18,3),'LineStyle','-','LineWidth',2,'Color',[.9,.9,0]);
hold on; plot(ControlSpectraAll(10:18,4),'LineStyle','-','LineWidth',2,'Color',[1,0,0]);
hold on; plot(ControlSpectraAll(10:18,5),'LineStyle','-','LineWidth',2,'Color',[1,0,1]);
hold on; plot(ControlSpectraAll(10:18,6),'LineStyle','-','LineWidth',2,'Color',[.5,.5,.5]);
hold on; plot(ControlSpectraAll(10:18,7),'LineStyle','-','LineWidth',2,'Color',[0,1,1]);
hold on; plot(ControlSpectraAll(10:18,8),'LineStyle','-','LineWidth',2,'Color',[0,0,0]);

hold on; plot(A0(10:18,1),'LineStyle','--','LineWidth',1,'Color',[0,0,1]);
hold on; plot(A0(10:18,2),'LineStyle','--','LineWidth',1,'Color',[0,1,0]);
hold on; plot(A0(10:18,3),'LineStyle','--','LineWidth',1,'Color',[.9,.9,0]);
hold on; plot(A0(10:18,4),'LineStyle','--','LineWidth',1,'Color',[1,0,0]);
hold on; plot(A0(10:18,5),'LineStyle','--','LineWidth',1,'Color',[1,0,1]);
hold on; plot(A0(10:18,6),'LineStyle','--','LineWidth',1,'Color',[.5,.5,.5]);
hold on; plot(A0(10:18,7),'LineStyle','--','LineWidth',1,'Color',[0,1,1]);
hold on; plot(A0(10:18,8),'LineStyle','--','LineWidth',1,'Color',[0,0,0]);

hold on; plot(AnmfRI(10:18,1),'LineStyle','-','Marker','*','LineWidth',1,'Color',[0,0,1]);
hold on; plot(AnmfRI(10:18,2),'LineStyle','-','Marker','*','LineWidth',1,'Color',[0,1,0]);
hold on; plot(AnmfRI(10:18,3),'LineStyle','-','Marker','*','LineWidth',1,'Color',[.9,.9,0]);
hold on; plot(AnmfRI(10:18,4),'LineStyle','-','Marker','*','LineWidth',1,'Color',[1,0,0]);
hold on; plot(AnmfRI(10:18,5),'LineStyle','-','Marker','*','LineWidth',1,'Color',[1,0,1]);
hold on; plot(AnmfRI(10:18,6),'LineStyle','-','Marker','*','LineWidth',1,'Color',[.5,.5,.5]);
hold on; plot(AnmfRI(10:18,7),'LineStyle','-','Marker','*','LineWidth',1,'Color',[0,1,1]);
hold on; plot(AnmfRI(10:18,8),'LineStyle','-','Marker','*','LineWidth',1,'Color',[0,0,0]);
set(ax,'XTick',1:7);
set(ax,'XTickLabel',[520,540,560,580,600,620,640,660,680])
xlabel('Center wavelength of spectral channels (nm)')
ylabel('Intensity (Normalized)')
axis tight;
box on;

%% figure 3
figure; title('Spectra on Cy3 Filter');  ax=gca;
hold on; plot(ControlSpectraAll(19:26,1),'LineStyle','-','LineWidth',2,'Color',[0,0,1]);
hold on; plot(ControlSpectraAll(19:26,2),'LineStyle','-','LineWidth',2,'Color',[0,1,0]);
hold on; plot(ControlSpectraAll(19:26,3),'LineStyle','-','LineWidth',2,'Color',[.9,.9,0]);
hold on; plot(ControlSpectraAll(19:26,4),'LineStyle','-','LineWidth',2,'Color',[1,0,0]);
hold on; plot(ControlSpectraAll(19:26,5),'LineStyle','-','LineWidth',2,'Color',[1,0,1]);
hold on; plot(ControlSpectraAll(19:26,6),'LineStyle','-','LineWidth',2,'Color',[.5,.5,.5]);
hold on; plot(ControlSpectraAll(19:26,7),'LineStyle','-','LineWidth',2,'Color',[0,1,1]);
hold on; plot(ControlSpectraAll(19:26,8),'LineStyle','-','LineWidth',2,'Color',[0,0,0]);

hold on; plot(A0(19:26,1),'LineStyle','--','LineWidth',1,'Color',[0,0,1]);
hold on; plot(A0(19:26,2),'LineStyle','--','LineWidth',1,'Color',[0,1,0]);
hold on; plot(A0(19:26,3),'LineStyle','--','LineWidth',1,'Color',[.9,.9,0]);
hold on; plot(A0(19:26,4),'LineStyle','--','LineWidth',1,'Color',[1,0,0]);
hold on; plot(A0(19:26,5),'LineStyle','--','LineWidth',1,'Color',[1,0,1]);
hold on; plot(A0(19:26,6),'LineStyle','--','LineWidth',1,'Color',[.5,.5,.5]);
hold on; plot(A0(19:26,7),'LineStyle','--','LineWidth',1,'Color',[0,1,1]);
hold on; plot(A0(19:26,8),'LineStyle','--','LineWidth',1,'Color',[0,0,0]);

hold on; plot(AnmfRI(19:26,1),'LineStyle','-','Marker','*','LineWidth',1,'Color',[0,0,1]);
hold on; plot(AnmfRI(19:26,2),'LineStyle','-','Marker','*','LineWidth',1,'Color',[0,1,0]);
hold on; plot(AnmfRI(19:26,3),'LineStyle','-','Marker','*','LineWidth',1,'Color',[.9,.9,0]);
hold on; plot(AnmfRI(19:26,4),'LineStyle','-','Marker','*','LineWidth',1,'Color',[1,0,0]);
hold on; plot(AnmfRI(19:26,5),'LineStyle','-','Marker','*','LineWidth',1,'Color',[1,0,1]);
hold on; plot(AnmfRI(19:26,6),'LineStyle','-','Marker','*','LineWidth',1,'Color',[.5,.5,.5]);
hold on; plot(AnmfRI(19:26,7),'LineStyle','-','Marker','*','LineWidth',1,'Color',[0,1,1]);
hold on; plot(AnmfRI(19:26,8),'LineStyle','-','Marker','*','LineWidth',1,'Color',[0,0,0]);
set(ax,'XTick',1:8);
set(ax,'XTickLabel',[570,590,610,630,650,670,690,710])
xlabel('Center wavelength of spectral channels (nm)')
ylabel('Intensity (Normalized)')
axis tight; 
box on;

%% figure 4
figure; title('Spectra on Texas Red and Cy5 Filter');  ax=gca;
hold on; plot(ControlSpectraAll(27:32,1),'LineStyle','-','LineWidth',2,'Color',[0,0,1]);
hold on; plot(ControlSpectraAll(27:32,2),'LineStyle','-','LineWidth',2,'Color',[0,1,0]);
hold on; plot(ControlSpectraAll(27:32,3),'LineStyle','-','LineWidth',2,'Color',[.9,.9,0]);
hold on; plot(ControlSpectraAll(27:32,4),'LineStyle','-','LineWidth',2,'Color',[1,0,0]);
hold on; plot(ControlSpectraAll(27:32,5),'LineStyle','-','LineWidth',2,'Color',[1,0,1]);
hold on; plot(ControlSpectraAll(27:32,6),'LineStyle','-','LineWidth',2,'Color',[.5,.5,.5]);
hold on; plot(ControlSpectraAll(27:32,7),'LineStyle','-','LineWidth',2,'Color',[0,1,1]);
hold on; plot(ControlSpectraAll(27:32,8),'LineStyle','-','LineWidth',2,'Color',[0,0,0]);

hold on; plot(A0(27:32,1),'LineStyle','--','LineWidth',1,'Color',[0,0,1]);
hold on; plot(A0(27:32,2),'LineStyle','--','LineWidth',1,'Color',[0,1,0]);
hold on; plot(A0(27:32,3),'LineStyle','--','LineWidth',1,'Color',[.9,.9,0]);
hold on; plot(A0(27:32,4),'LineStyle','--','LineWidth',1,'Color',[1,0,0]);
hold on; plot(A0(27:32,5),'LineStyle','--','LineWidth',1,'Color',[1,0,1]);
hold on; plot(A0(27:32,6),'LineStyle','--','LineWidth',1,'Color',[.5,.5,.5]);
hold on; plot(A0(27:32,7),'LineStyle','--','LineWidth',1,'Color',[0,1,1]);
hold on; plot(A0(27:32,8),'LineStyle','--','LineWidth',1,'Color',[0,0,0]);

hold on; plot(AnmfRI(27:32,1),'LineStyle','-','Marker','*','LineWidth',1,'Color',[0,0,1]);
hold on; plot(AnmfRI(27:32,2),'LineStyle','-','Marker','*','LineWidth',1,'Color',[0,1,0]);
hold on; plot(AnmfRI(27:32,3),'LineStyle','-','Marker','*','LineWidth',1,'Color',[.9,.9,0]);
hold on; plot(AnmfRI(27:32,4),'LineStyle','-','Marker','*','LineWidth',1,'Color',[1,0,0]);
hold on; plot(AnmfRI(27:32,5),'LineStyle','-','Marker','*','LineWidth',1,'Color',[1,0,1]);
hold on; plot(AnmfRI(27:32,6),'LineStyle','-','Marker','*','LineWidth',1,'Color',[.5,.5,.5]);
hold on; plot(AnmfRI(27:32,7),'LineStyle','-','Marker','*','LineWidth',1,'Color',[0,1,1]);
hold on; plot(AnmfRI(27:32,8),'LineStyle','-','Marker','*','LineWidth',1,'Color',[0,0,0]);
set(ax,'XTick',1:6);
set(ax,'XTickLabel',[580,590,620,640,660,680])
xlabel('Center wavelength of spectral channels (nm)')
ylabel('Intensity (Normalized)')
axis tight;
box on;

%% figure 5
figure; ax=gca;
hold on; plot(ControlSpectraAll(33:end,1),'LineStyle','-','LineWidth',2,'Color',[0,0,1]);
hold on; plot(ControlSpectraAll(33:end,2),'LineStyle','-','LineWidth',2,'Color',[0,1,0]);
hold on; plot(ControlSpectraAll(33:end,3),'LineStyle','-','LineWidth',2,'Color',[.9,.9,0]);
hold on; plot(ControlSpectraAll(33:end,4),'LineStyle','-','LineWidth',2,'Color',[1,0,0]);
hold on; plot(ControlSpectraAll(33:end,5),'LineStyle','-','LineWidth',2,'Color',[1,0,1]);
hold on; plot(ControlSpectraAll(33:end,6),'LineStyle','-','LineWidth',2,'Color',[.5,.5,.5]);
hold on; plot(ControlSpectraAll(33:end,7),'LineStyle','-','LineWidth',2,'Color',[0,1,1]);
hold on; plot(ControlSpectraAll(33:end,8),'LineStyle','-','LineWidth',2,'Color',[0,0,0]);

hold on; plot(A0(33:end,1),'LineStyle','--','LineWidth',1,'Color',[0,0,1]);
hold on; plot(A0(33:end,2),'LineStyle','--','LineWidth',1,'Color',[0,1,0]);
hold on; plot(A0(33:end,3),'LineStyle','--','LineWidth',1,'Color',[.9,.9,0]);
hold on; plot(A0(33:end,4),'LineStyle','--','LineWidth',1,'Color',[1,0,0]);
hold on; plot(A0(33:end,5),'LineStyle','--','LineWidth',1,'Color',[1,0,1]);
hold on; plot(A0(33:end,6),'LineStyle','--','LineWidth',1,'Color',[.5,.5,.5]);
hold on; plot(A0(33:end,7),'LineStyle','--','LineWidth',1,'Color',[0,1,1]);
hold on; plot(A0(33:end,8),'LineStyle','--','LineWidth',1,'Color',[0,0,0]);

hold on; plot(AnmfRI(33:end,1),'LineStyle','-','Marker','*','LineWidth',1,'Color',[0,0,1]);
hold on; plot(AnmfRI(33:end,2),'LineStyle','-','Marker','*','LineWidth',1,'Color',[0,1,0]);
hold on; plot(AnmfRI(33:end,3),'LineStyle','-','Marker','*','LineWidth',1,'Color',[.9,.9,0]);
hold on; plot(AnmfRI(33:end,4),'LineStyle','-','Marker','*','LineWidth',1,'Color',[1,0,0]);
hold on; plot(AnmfRI(33:end,5),'LineStyle','-','Marker','*','LineWidth',1,'Color',[1,0,1]);
hold on; plot(AnmfRI(33:end,6),'LineStyle','-','Marker','*','LineWidth',1,'Color',[.5,.5,.5]);
hold on; plot(AnmfRI(33:end,7),'LineStyle','-','Marker','*','LineWidth',1,'Color',[0,1,1]);
hold on; plot(AnmfRI(33:end,8),'LineStyle','-','Marker','*','LineWidth',1,'Color',[0,0,0]);
set(ax,'XTick',1:3);
set(ax,'XTickLabel',[680,700,720])
xlabel('Center wavelength of spectral channels (nm)')
ylabel('Intensity (Normalized)')
axis tight;
box on;


end