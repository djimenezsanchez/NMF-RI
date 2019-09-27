function [] = displayResult(AnmfRI,ControlSpectraAll,A0)

%% figure 1
figure;  ax=gca;
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

hold on; plot(AnmfRI(1:9,1),'LineStyle','-','Marker','*','MarkerSize',9,'LineWidth',1,'Color',[0,0,1]);
hold on; plot(AnmfRI(1:9,2),'LineStyle','-','Marker','*','MarkerSize',9,'LineWidth',1,'Color',[0,1,0]);
hold on; plot(AnmfRI(1:9,3),'LineStyle','-','Marker','*','MarkerSize',9,'LineWidth',1,'Color',[.9,.9,0]);
hold on; plot(AnmfRI(1:9,4),'LineStyle','-','Marker','*','MarkerSize',9,'LineWidth',1,'Color',[1,0,0]);
hold on; plot(AnmfRI(1:9,5),'LineStyle','-','Marker','*','MarkerSize',9,'LineWidth',1,'Color',[1,0,1]);
hold on; plot(AnmfRI(1:9,6),'LineStyle','-','Marker','*','MarkerSize',9,'LineWidth',1,'Color',[.5,.5,.5]);
hold on; plot(AnmfRI(1:9,7),'LineStyle','-','Marker','*','MarkerSize',9,'LineWidth',1,'Color',[0,1,1]);
hold on; plot(AnmfRI(1:9,8),'LineStyle','-','Marker','*','MarkerSize',9,'LineWidth',1,'Color',[0,0,0]);
% 
% hold on; plot(Anmf_withoutAF(1:9,1),'LineStyle','-','Marker','*','MarkerSize',9,'LineWidth',1,'Color',[0,0,1]);
% hold on; plot(Anmf_withoutAF(1:9,2),'LineStyle','-','Marker','*','MarkerSize',9,'LineWidth',1,'Color',[0,1,0]);
% % hold on; plot(Anmf_withoutAF(1:9,3),'LineStyle','-','Marker','*','MarkerSize',9,'LineWidth',1,'Color',[.9,.9,0]);
% % hold on; plot(Anmf_withoutAF(1:9,4),'LineStyle','-','Marker','*','MarkerSize',9,'LineWidth',1,'Color',[1,0,0]);
% % hold on; plot(Anmf_withoutAF(1:9,5),'LineStyle','-','Marker','*','MarkerSize',9,'LineWidth',1,'Color',[1,0,1]);
% % hold on; plot(Anmf_withoutAF(1:9,6),'LineStyle','-','Marker','*','MarkerSize',9,'LineWidth',1,'Color',[.5,.5,.5]);
% % hold on; plot(Anmf_withoutAF(1:9,7),'LineStyle','-','Marker','*','MarkerSize',9,'LineWidth',1,'Color',[0,1,1]);
% hold on; plot(zeros(1,9),'LineStyle','-','Marker','*','MarkerSize',9,'LineWidth',1,'Color',[0,0,0]);
% 
% hold on; plot(Anmf_withstraightAF(1:9,1),'LineStyle','-','Marker','+','MarkerSize',9,'LineWidth',1,'Color',[0,0,1]);
% hold on; plot(Anmf_withstraightAF(1:9,2),'LineStyle','-','Marker','+','MarkerSize',9,'LineWidth',1,'Color',[0,1,0]);
% % hold on; plot(Anmf_withstraightAF(1:9,3),'LineStyle','-','Marker','+','MarkerSize',9,'LineWidth',1,'Color',[.9,.9,0]);
% % hold on; plot(Anmf_withstraightAF(1:9,4),'LineStyle','-','Marker','+','MarkerSize',9,'LineWidth',1,'Color',[1,0,0]);
% % hold on; plot(Anmf_withstraightAF(1:9,5),'LineStyle','-','Marker','+','MarkerSize',9,'LineWidth',1,'Color',[1,0,1]);
% % hold on; plot(Anmf_withstraightAF(1:9,6),'LineStyle','-','Marker','+','MarkerSize',9,'LineWidth',1,'Color',[.5,.5,.5]);
% % hold on; plot(Anmf_withstraightAF(1:9,7),'LineStyle','-','Marker','+','MarkerSize',9,'LineWidth',1,'Color',[0,1,1]);
% hold on; plot(Anmf_withstraightAF(1:9,8),'LineStyle','-','Marker','+','MarkerSize',9,'LineWidth',1,'Color',[0,0,0]);

h = legend('DAPI Reference', 'Opal 520 Reference', 'Opal 540 Reference', 'Opal 570 Reference', ...
    'Opal 620 Reference','Opal 650 Reference', 'Opal 690 Reference', 'AF Reference', ...
    'DAPI Initialization', 'Opal 520 Initialization', 'Opal 540 Initialization', 'Opal 570 Initialization', ...
    'Opal 620 Initialization','Opal 650 Initialization', 'Opal 690 Initialization', 'AF Initialization', ...
    'DAPI NMF-RI', 'Opal 520 NMF-RI', 'Opal 540 NMF-RI', 'Opal 570 NMF-RI', ...
    'Opal 620 NMF-RI','Opal 650 NMF-RI', 'Opal 690 NMF-RI', 'AF NMF-RI' ...
    );
%     'DAPI NMF-RI (Without AF)', 'Opal 520 NMF-RI (Without AF)',...
%     'Opal 540 NMF-RI (Without AF)', 'Opal 570 NMF-RI (Without AF)', ...
%     'Opal 620 NMF-RI (Without AF)','Opal 650 NMF-RI (Without AF)',...
%     'Opal 690 NMF-RI (Without AF)', 'AF NMF-RI (Without AF)', ...
%     'DAPI NMF-RI (Without AF reference)', 'Opal 520 NMF-RI (Without AF reference)',...
%     'Opal 540 NMF-RI (Without AF reference)', 'Opal 570 NMF-RI (Without AF reference)', ...
%     'Opal 620 NMF-RI (Without AF reference)','Opal 650 NMF-RI (Without AF reference)',...
%     'Opal 690 NMF-RI (Without AF reference)', 'AF NMF-RI (Without AF reference)');
h.NumColumns = 3;
% hold on; plot(Anmf_withoutAF(1:9,8),'LineStyle','-','Marker','o','LineWidth',1,'Color',[0,0,0]);
set(ax,'XTick',1:9);
set(ax,'XTickLabel',[440,460,480,500,520,540,560,580,600])
xlabel('Center wavelength of spectral channels (nm)')
ylabel('Intensity (Normalized)')
axis tight;
box on;

%% figure 2
figure;  ax=gca;
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

hold on; plot(AnmfRI(10:18,1),'LineStyle','-','Marker','*','MarkerSize',9,'LineWidth',1,'Color',[0,0,1]);
hold on; plot(AnmfRI(10:18,2),'LineStyle','-','Marker','*','MarkerSize',9,'LineWidth',1,'Color',[0,1,0]);
hold on; plot(AnmfRI(10:18,3),'LineStyle','-','Marker','*','MarkerSize',9,'LineWidth',1,'Color',[.9,.9,0]);
hold on; plot(AnmfRI(10:18,4),'LineStyle','-','Marker','*','MarkerSize',9,'LineWidth',1,'Color',[1,0,0]);
hold on; plot(AnmfRI(10:18,5),'LineStyle','-','Marker','*','MarkerSize',9,'LineWidth',1,'Color',[1,0,1]);
hold on; plot(AnmfRI(10:18,6),'LineStyle','-','Marker','*','MarkerSize',9,'LineWidth',1,'Color',[.5,.5,.5]);
hold on; plot(AnmfRI(10:18,7),'LineStyle','-','Marker','*','MarkerSize',9,'LineWidth',1,'Color',[0,1,1]);
hold on; plot(AnmfRI(10:18,8),'LineStyle','-','Marker','*','MarkerSize',9,'LineWidth',1,'Color',[0,0,0]);

% hold on; plot(Anmf_withoutAF(10:18,1),'LineStyle','-','Marker','*','MarkerSize',9,'LineWidth',1,'Color',[0,0,1]);
% hold on; plot(Anmf_withoutAF(10:18,2),'LineStyle','-','Marker','*','MarkerSize',9,'LineWidth',1,'Color',[0,1,0]);
% % hold on; plot(Anmf_withoutAF(10:18,3),'LineStyle','-','Marker','*','MarkerSize',9,'LineWidth',1,'Color',[.9,.9,0]);
% % hold on; plot(Anmf_withoutAF(10:18,4),'LineStyle','-','Marker','*','MarkerSize',9,'LineWidth',1,'Color',[1,0,0]);
% % hold on; plot(Anmf_withoutAF(10:18,5),'LineStyle','-','Marker','*','MarkerSize',9,'LineWidth',1,'Color',[1,0,1]);
% % hold on; plot(Anmf_withoutAF(10:18,6),'LineStyle','-','Marker','*','MarkerSize',9,'LineWidth',1,'Color',[.5,.5,.5]);
% % hold on; plot(Anmf_withoutAF(10:18,7),'LineStyle','-','Marker','*','MarkerSize',9,'LineWidth',1,'Color',[0,1,1]);
% hold on; plot(zeros(9,1),'LineStyle','-','Marker','*','LineWidth',1,'Color',[0,0,0]);
% 
% hold on; plot(Anmf_withstraightAF(10:18,1),'LineStyle','-','Marker','+','MarkerSize',9,'LineWidth',1,'Color',[0,0,1]);
% hold on; plot(Anmf_withstraightAF(10:18,2),'LineStyle','-','Marker','+','MarkerSize',9,'LineWidth',1,'Color',[0,1,0]);
% % hold on; plot(Anmf_withstraightAF(10:18,3),'LineStyle','-','Marker','+','MarkerSize',9,'LineWidth',1,'Color',[.9,.9,0]);
% % hold on; plot(Anmf_withstraightAF(10:18,4),'LineStyle','-','Marker','+','MarkerSize',9,'LineWidth',1,'Color',[1,0,0]);
% % hold on; plot(Anmf_withstraightAF(10:18,5),'LineStyle','-','Marker','+','MarkerSize',9,'LineWidth',1,'Color',[1,0,1]);
% % hold on; plot(Anmf_withstraightAF(10:18,6),'LineStyle','-','Marker','+','MarkerSize',9,'LineWidth',1,'Color',[.5,.5,.5]);
% % hold on; plot(Anmf_withstraightAF(10:18,7),'LineStyle','-','Marker','+','MarkerSize',9,'LineWidth',1,'Color',[0,1,1]);
% hold on; plot(Anmf_withstraightAF(10:18,8),'LineStyle','-','Marker','+','MarkerSize',9,'LineWidth',1,'Color',[0,0,0]);

set(ax,'XTick',1:7);
set(ax,'XTickLabel',[520,540,560,580,600,620,640,660,680])
xlabel('Center wavelength of spectral channels (nm)')
ylabel('Intensity (Normalized)')
axis tight;
box on;

end