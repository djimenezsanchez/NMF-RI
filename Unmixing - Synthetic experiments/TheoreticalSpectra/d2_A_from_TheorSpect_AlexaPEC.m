%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain ncx2 crosstalk matrix from Alexa488 and RPE theoretical spectra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = d2_A_from_TheorSpect_RPEPEC(nc)

% Datos de interes
 warning('off', 'MATLAB:interp1:UsePCHIP')
% % % canales = nc;
% % % posicioninicio = 17; %en el vector cwl (este dato depende del inicio de los fluorocromos a estudair)
imFolder = 'TheorSpect/';
    

% Cargamos ahora las curva teóricas y las interpolamos en todo el rango

ALEXA_T=load([imFolder 'alexa488.txt'])';
ALEXA_T_X=min(ALEXA_T(1,:)):1:max(ALEXA_T(1,:));
ALEXA_T_interp=interp1(ALEXA_T(1,:),ALEXA_T(2,:),ALEXA_T_X,'cubic');

PEC_T=load([imFolder 'PECyN5.txt'])';
PEC_T(2,:)=PEC_T(2,:)/max(PEC_T(2,:));
PEC_T_X=min(PEC_T(1,:)):1:max(PEC_T(1,:));
PEC_T_interp=interp1(PEC_T(1,:),PEC_T(2,:),PEC_T_X,'cubic');


minX = min([ALEXA_T_X, PEC_T_X]);
maxX = max([ALEXA_T_X, PEC_T_X]);
minX = 460;
maxX = 780;
AllX = minX:0.01:maxX;

ALEXA_all = interp1(ALEXA_T(1,:),ALEXA_T(2,:),AllX,'cubic',0);
PEC_all = interp1(PEC_T(1,:),PEC_T(2,:),AllX,'cubic',0);

step = floor(length(AllX)/nc);
% step = (maxX-minX)/nc;
% temp = minX;
for i=1:nc
%     pos(i) = mean([temp temp+step]);
%     temp = temp+step;
    A(i,1) = sum( ALEXA_all( (step*(i-1)+1):step*i ) );
    A(i,2) = sum( PEC_all( (step*(i-1)+1):step*i ) );
end

% 
% figure, plot(AllX, ALEXA_all, 'g'),
% hold on, plot(AllX, PEC_all, 'r'),
% plot([AllX(1) AllX(1)], [0 1], 'b')
% for i=1:nc
%     plot([AllX(i*step) AllX(i*step)], [0 1], 'b')
% end
% title('Theoretical Spectra');
% legend('ALEXA 488', 'RPE');




end
