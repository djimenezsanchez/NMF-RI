%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain ncx2 crosstalk matrix from Alexa488 and RPE theoretical spectra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = d2_A_from_TheorSpect_AlexaRPEPEC(nc)

% Datos de interes
 warning('off', 'MATLAB:interp1:UsePCHIP')
% % % canales = nc;
% % % posicioninicio = 17; %en el vector cwl (este dato depende del inicio de los fluorocromos a estudair)
imFolder = 'TheorSpect/';
    

% Cargamos ahora las curva teóricas y las interpolamos en todo el rango
Alexa_T=load([imFolder 'alexa488.txt'])';
Alexa_T(2,:)=Alexa_T(2,:)/max(Alexa_T(2,:));
Alexa_T_X=min(Alexa_T(1,:)):1:max(Alexa_T(1,:));
Alexa_T_interp=interp1(Alexa_T(1,:),Alexa_T(2,:),Alexa_T_X,'cubic');

PEC_T=load([imFolder 'PECyN5.txt'])';
PEC_T_X=min(PEC_T(1,:)):1:max(PEC_T(1,:));
PEC_T_interp=interp1(PEC_T(1,:),PEC_T(2,:),PEC_T_X,'cubic');

RPE_T=load([imFolder 'RPE.txt'])';
RPE_T(2,:)=RPE_T(2,:)/max(RPE_T(2,:));
RPE_T_X=min(RPE_T(1,:)):1:max(RPE_T(1,:));
RPE_T_interp=interp1(RPE_T(1,:),RPE_T(2,:),RPE_T_X,'cubic');


minX = min([PEC_T_X, RPE_T_X, Alexa_T_X]);
maxX = max([PEC_T_X, RPE_T_X, Alexa_T_X]);
minX = 460;
maxX = 740;
AllX = minX:0.01:maxX;

PEC_all = interp1(PEC_T(1,:),PEC_T(2,:),AllX,'cubic',0);
RPE_all = interp1(RPE_T(1,:),RPE_T(2,:),AllX,'cubic',0);
Alexa_all = interp1(Alexa_T(1,:),Alexa_T(2,:),AllX,'cubic',0);
PEC_all = PEC_all./max(PEC_all);
RPE_all = RPE_all./max(RPE_all);
Alexa_all = Alexa_all./max(Alexa_all);

step = floor(length(AllX)/nc);
% step = (maxX-minX)/nc;
% temp = minX;
for i=1:nc
    A(i,3) = sum( PEC_all( (step*(i-1)+1):step*i ) );
    A(i,2) = sum( RPE_all( (step*(i-1)+1):step*i ) );
    A(i,1) = sum( Alexa_all( (step*(i-1)+1):step*i ) );
end


% figure, plot(AllX, PEC_all, 'g'),
% hold on, plot(AllX, RPE_all, 'r'), plot(AllX, Alexa_all, 'y')
% plot([AllX(1) AllX(1)], [0 1], 'b')
% for i=1:nc
%     plot([AllX(i*step) AllX(i*step)], [0 1], 'b')
% end
% title('Theoretical Spectra');
% legend('PEC', 'RPE', 'Alexa');
% 
% 
% 

end
