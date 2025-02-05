%% ----------------- 1) CARGAR ARCHIVOS Y JUNTARLOS
% Concatena archivos abf
clear; clc; close all
tic
abfsCell = dir('*.abf');
list_abf  = [];
for nabf = 1:length(abfsCell)
    list_abf = [list_abf;abfsCell(nabf).name];
end
[data,fs,files_spec] = abf_cat(list_abf(1:end,:)); % concatenate channel throught all files
%[data,fs,files_spec] = abf_cat(list_abf(3:3,:)); % concatenate channel throught all files


%abfsCell = dir('*.abf');
abfsCell = ls('*.abf');
[data,fs,files_spec] = abf_cat(abfsCell(1:end,:)); % concatenate channel throught all files

%% Quitar ruido 50 Hz (notch filter)
Fs= 20000;
wo = 50/(Fs/2);   factor  = 10;
bw= wo/ factor;
[b, a] = iirnotch(wo, bw);
data_filt = filtfilt(b,a,data);    data= (data_filt);
clear data_filt
%% ----------------- 2) REDUCIR FRECUENCIA MUESTREO (2000 Hz)
% Esto permite ejecutar mas rapido el analisis 
% Para el analisis de ripples usar la señal original 
factor_ds = 10;
disp('downsampling...'); data_ds = downsample(data,factor_ds);
fs_ds = Fs/factor_ds;

%% ----------------- 3) GRABAR
disp('saving...'); save('raw_WF.mat','data','fs','data_ds','fs_ds','-v7.3');
%save('raw_WF.mat','data_ds','fs_ds','-append');
% clear data fs
toc

%% ----------------- 4) GENERAR ESPECTRO EN CADA CANAL
% Espectro de potencia para diferentes bandas
% viewLFP 3 argumentos: matriz a dibujar, tamaño de ventana (en seg),
% expansion del canal
viewLFP(data_ds(:,4:16),fs_ds,6,3); % Visualizar registro

% Espectro total. Pregunta por el canal de max de theta (por si midiera algo raro)
% Decir y
% Si quito el ,1) del argumento final, no me pregunta por el max theta ni saca timespec
% Pongo data_ds(:, 4:16) porque los tres primeros canales no tienen LFP
% Por tanto, a las cifras de numero de canal del grafico hay que sumar 3
[pS,band,freq]=powerSpectrumProfileAlb(data_ds(:,1:16),fs_ds,1); 
[pS,band,freq]=powerSpectrumProfileAlb_noSD(data_ds(:,1:16),fs_ds,1); 

ylim([40,75])
%% ----------------- 5) SELECCIONAR CANALES Y ZONAS ESPECIALES
ch.pyr = 11; % Suele ser el 1er canal en que aparece max local en HFOs y min local en Theta
ch.rad = [13]; % Intermedio entre pyr y slm
ch.slm = [14]; % Suele coincidir con el max de theta, continuando la curva desde pyr
ch.hil = []; % Lo busco a ojo por las espigas de dentado 
ch.wheel = []; % No se grabo. Si grabara movimiento de rueda


% Para estos registros, dejar todos los diary vacios (como mucho alguna zona en discard)
diary.dmso = []; % in seconds. Zonas en las que se añadio DMSO
diary.PTZ_CNO = []; % in seconds. Zonas en las que se añadio CNO
diary.EP = []; % in seconds. Zonas en las que se dieron pulsos de electroporacion
diary.discard = []; % in seconds. Zonas a descartar

%% ---------------- 6) GRABAR
save('raw_WF.mat','ch','diary','band','freq','-append');

%% ---------------- 7)  ELIMINAR ZONAS DE DESCARTE (NO USAR)
% data_c = dataCleaner(data,fs,diary); % TAREA 1: QUIZÁ HAY QUE CAMBIAR LA DURACIÓN DE LA VENTANA QUE BLANQUEAMOS!!!
% data_ds_c = dataCleaner(data_ds,fs_ds,diary);

%% ---------------- 8) CALCULAR ZONAS EN THETA Y LIA (large-amplitude irregular activity)
% Basicamente theta es cuando corre y LIA cuando no corre (es aqui cuando aparecen ripples)

% Ratios Theta/Delta (alto cuando corre) y Theta/Beta (bajo cuando corre) y
% tiempos de cada ratio --> Se usara solo el theta/delta
time_ds = linspace(0,size(data_ds, 1)/fs_ds, size(data_ds, 1))';
[thetaDeltaStr thetaBetaStr]= thetaFinderWheel(data_ds(:,ch.slm), fs_ds, 2, time_ds); % find theta instants, change theta/delta threslhold to refine

% Calcula el espectro de potencias general
[S, S_det, t, f] = powerSpectrumWide(data_ds, fs_ds, thetaDeltaStr, ch);

freq_plot = find(f<30);  % Pone limite a la frecuencia a representar
fs_spc = (length(data_ds)/fs_ds)/length(thetaDeltaStr.seg); 
theta_idx = find(thetaDeltaStr.seg);
t_theta = length(theta_idx)*fs_spc;
LIA_idx = find(~thetaDeltaStr.seg);
t_LIA = length(LIA_idx)*fs_spc;

% powerSpectrum THETA vs. LIA: theta & gammaS/gammaF/gammaSuperF
bandCon = powerSpecCh(fs_ds, data_ds, theta_idx, LIA_idx);

% Plots de los espectros de frecuencia
figPowerSpec(t, f, S_det, ch, [-15, 15], freq_plot, 1:length(thetaDeltaStr.seg), 'General');
figPowerSpec([1 t_theta], f, S_det, ch, [-15, 15], freq_plot, theta_idx, 'Theta');
figPowerSpec([1 t_LIA], f, S_det, ch, [-15, 15], freq_plot, LIA_idx, 'LIA');

figPowerSpecFreq(f, S, ch, theta_idx, LIA_idx, 'Theta', 'LIA');
figPowerSpecFreq(f, S_det, ch, theta_idx, LIA_idx, 'Theta detrended', 'LIA detrended');

figure;
plot(bandCon.ch,zscore(bandCon.thetaCon.theta));

%% ---------------- 9) GRABAR
save('raw_WF.mat','thetaDeltaStr','thetaBetaStr','S','S_det','t','f', 'bandCon', '-append');
save('raw_WF.mat','thetaDeltaStr','thetaBetaStr', '-append');


%% ------------------ 10) Find ripples
% Añade algunos canales adicionales para la funcion de buscar ripples
% data =data(1:1200*fs,:);
ch.CA1=ch.pyr;
ch.noise=3;
ch.CA3=[];
ch.opto=[];
% Argumentos: matriz de datos, canales, canal de piramidal, freq muestreo,
% intervalo de frecuencias a detectar, umbral de deteccion, intervalo para
% calcular media y desviacion de la envolvente de señal filtrada (sobre ella se buscan los ripples en funcion de un umbral)
data_filt(770*fs:end, ch.pyr) = data_filt(770*fs:end, ch.pyr-1);
%[spontaneous_ripples, media, SD] = rippleFinderAttila(data, ch, ch.pyr, fs, [90,250], 3,[1,45]); % 1,34  [9, [1, 2]]   % pruebo con 3
[spontaneous_ripples, media, SD] = rippleFinderNew(data, ch, 'chrip',ch.CA1,'fs',fs,'Fpass',[90, 250],'threshold',3,'SD', 0.005); % 'SD', 0.0024);
[media, SD] = mediaSDrippleFinder(data, ch, ch.pyr, fs, [90,250], 3);

% Propiedades de ripples
% añade cunciones a la estructura que son necesarias para la funcion
% propSpVerde --> calcula propiedades de los ripples
cell_name= 'Nuria54_C_1'; % 'Nuria54_E_2'  Thy1xGCaMP_4_2i
spontaneous_ripples.media= ones(size(spontaneous_ripples.mat_CA1,2),1).*media;
spontaneous_ripples.estandar= ones(size(spontaneous_ripples.mat_CA1,2),1).*SD;
spontaneous_ripples.filename=cell(size(spontaneous_ripples.mat_CA1,2),1);
spontaneous_ripples.num_evento = [];
for i=1:size(spontaneous_ripples.mat_CA1,2)
    spontaneous_ripples.filename{i,1}= cell_name;
    spontaneous_ripples.num_evento(i,1)= i;
end
spontaneous_ripples= propSpVerde(spontaneous_ripples, fs);
spontaneous_ripples.frequency= spontaneous_ripples.valores{:,5};
spontaneous_ripples.powerSD= spontaneous_ripples.valores{:,3};
%

rip_ts=spontaneous_ripples.time_ripples;
mat_rip=spontaneous_ripples.mat_CA1;

[ripSpec]=ripSpectrogramManu2(mat_rip,fs); ripSpec.rip_ts = rip_ts; 

save('raw_WF.mat','ripSpec', 'ch', 'spontaneous_ripples','-append'); %% SAVE POINT!!


%% Quitar ripples
media = spontaneous_ripples.media(1);
SD = spontaneous_ripples.estandar(1);
spontaneous_ripples = quitaripples(spontaneous_ripples, fs);

%%
% Plot ripples
%ripPlot;
th_theta = 3; th_LIA = 1.8; % th_theta = 3; th_LIA = 1.8;

[pStemp.nodB, pStemp.dB] = pStemporalWinMat(data_ds, fs_ds, ch.pyr-1, [], 2, [1, 1000], true);
time2 = linspace(0, size(data_ds,1)/fs_ds -1.5, size(pStemp.nodB, 2))';
ratio = nanmean(pStemp.nodB(7:19,:)) ./ nanmean(pStemp.nodB(1:6,:)); % theta 5-8 / delta
ratioEnv = movmean( movmean( sgolayfilt(ratio,50,51), 20), 20); 
ratioEnv =(ratioEnv.^2)/2; %./2
figure; plot(time2, ratio); hold on; plot(time2, ratioEnv); yline(th_LIA,'r'); yline(th_theta,'g');
t_LIACo = length(ratioEnv(ratioEnv < th_LIA))./length(ratioEnv).*time2(end); 
t_ThetaCo = length(ratioEnv(ratioEnv > th_theta))./length(ratioEnv).*time2(end); 

ripCo = spontaneous_ripples;
abundCo= length(ripCo.idx_ripples);
abundCo_lento= length(ripCo.valores{ripCo.valores.Frequencia<120,1});
abundCoHz= abundCo/t_LIACo;
abundCoHz_lento= abundCo_lento/t_LIACo;
% ripPlotCSD;

%% upsample
ratiolong = nan(size(data,1),1);
timelong = round(linspace(0,size(data, 1)/fs, size(data, 1))', 4);
timeshort = round(time2,1);

for i = 1: length(time2)-1
    if any([10:10:100] == round(i / (length(time2)-1)*100,0))
    	fprintf('%d \n', round(i / (length(time2)-1)*100,0));
    end
    a1 = find(timelong == timeshort(i), 1, 'first');
    a2 = find(timelong == timeshort(i+1) - 0.0001, 1, 'last');
    ratiolong(a1: a2) = ratioEnv(i);
end
ratiolong(ratiolong > 2) = nan;  %
data2 = data;
data2(isnan(ratiolong), :) = 0;
[spontaneous_ripples, media, SD] = rippleFinderNew(data2, ch, 'chrip',ch.CA1,'fs',fs,'Fpass',[90, 250],'threshold',3,'SD', SD);

% Abundancia de ripples
if isempty(ch.slm)
    ch.slm = 16;
end

% data_ds =data_ds(1:1200*fs_ds,:);
time_ds = linspace(0,size(data_ds, 1)/fs_ds, size(data_ds, 1))';

[thetaDeltaStr thetaBetaStr]= thetaFinderWheel(data_ds(:,ch.slm), fs_ds, 2, time_ds); % find theta instants, change theta/delta threslhold to refine

sectoresCo= thetaDeltaStr.seg;
fsCo = size(data_ds,1)/length(sectoresCo);
t_thetaCo = length(find(sectoresCo))*0.25;
t_LIACo = length(find(~sectoresCo))*0.25;
ripCo = spontaneous_ripples;
abundCo= length(ripCo.idx_ripples);
abundCo_lento= length(ripCo.valores{ripCo.valores.Frequencia<120,1});
abundCoHz= abundCo/t_LIACo;
abundCoHz_lento= abundCo_lento/t_LIACo;
%




%% Buscar ciclos de theta
% Usamos el downsampled para buscar
chTheta = ch.slm;
[iCyc, iMiddle] = KA_find_thetaAlb(data_ds, fs_ds, [], ch, chTheta);
[thetaSpec] = ThetaCyclesLFP(iCyc,data_ds(:,ch.pyr-1),data_ds(:,ch.slm), fs_ds, 'none');
iCyc = iCyc * factor_ds;   iMiddle = iMiddle * factor_ds;
thetaSpec.all_cyc = thetaSpec.all_cyc * 10;
save('raw_WF.mat','iCyc', 'iMiddle', 'thetaSpec','-append'); %% SAVE POINT!!

%% ---------------- 11) Crear base de datos
% Transformar Inf en Nan en S
campo = fields (S);
for i=1 : length(campo)
    paso = S.(campo{i});
    paso(:, isinf(sum(paso))) = nan;
    S.(campo{i}) = paso; 
end
%

% Animal : Control / Tratado
ruta = pwd;
cd 'E:\Users\database'
load ovario.mat

animal = 'A7_C_4';  % Cambiar de 'Control' a 'Tratado'
% grupo = 'Tratado'; 
% grupo = 'Control';  / Tratado / Control2
% grupo = 'Control2';

spontaneous_ripples= rmfield(spontaneous_ripples,{'mat_light', 'mat_time', 'IFA'});
theta_idx = find(thetaDeltaStr.seg);
LIA_idx = find(~thetaDeltaStr.seg);

Ovx.(grupo).(animal).name = animal;
Ovx.(grupo).(animal).ruta = ruta;
Ovx.(grupo).(animal).Nraton = 2;
Ovx.(grupo).(animal).rip = spontaneous_ripples;
Ovx.(grupo).(animal).rip.ripSpec = ripSpec;
Ovx.(grupo).(animal).theta.thetaSpec = thetaSpec;
Ovx.(grupo).(animal).theta.iMiddle = iMiddle;
Ovx.(grupo).(animal).theta.iCyc = iCyc;
Ovx.(grupo).(animal).fs= fs;
Ovx.(grupo).(animal).fs_ds= fs_ds;
Ovx.(grupo).(animal).thetaBeta= thetaBetaStr;
Ovx.(grupo).(animal).thetaDelta= thetaDeltaStr;
Ovx.(grupo).(animal).ch= ch;
Ovx.(grupo).(animal).band.band= band;
Ovx.(grupo).(animal).band.bandCon= bandCon;

Ovx.(grupo).(animal).Spec.S_LIA.or.med= nanmean(S.or(:, LIA_idx),2);
Ovx.(grupo).(animal).Spec.S_LIA.pyr.med= nanmean(S.pyr(:, LIA_idx),2);
Ovx.(grupo).(animal).Spec.S_LIA.rad.med= nanmean(S.rad(:, LIA_idx),2);
Ovx.(grupo).(animal).Spec.S_LIA.slm.med= nanmean(S.slm(:, LIA_idx),2);
Ovx.(grupo).(animal).Spec.S_LIA.hil.med= nanmean(S.hil(:, LIA_idx),2);

Ovx.(grupo).(animal).Spec.S_theta.or.med= nanmean(S.or(:, theta_idx),2);
Ovx.(grupo).(animal).Spec.S_theta.pyr.med= nanmean(S.pyr(:, theta_idx),2);
Ovx.(grupo).(animal).Spec.S_theta.rad.med= nanmean(S.rad(:, theta_idx),2);
Ovx.(grupo).(animal).Spec.S_theta.slm.med= nanmean(S.slm(:, theta_idx),2);
Ovx.(grupo).(animal).Spec.S_theta.hil.med= nanmean(S.hil(:, theta_idx),2);

Ovx.(grupo).(animal).Spec.S_LIA.or.SD= nanstd(S.or(:, LIA_idx),[],2);
Ovx.(grupo).(animal).Spec.S_LIA.pyr.SD= nanstd(S.pyr(:, LIA_idx),[],2);
Ovx.(grupo).(animal).Spec.S_LIA.rad.SD= nanstd(S.rad(:, LIA_idx),[],2);
Ovx.(grupo).(animal).Spec.S_LIA.slm.SD= nanstd(S.slm(:, LIA_idx),[],2);
Ovx.(grupo).(animal).Spec.S_LIA.hil.SD= nanstd(S.hil(:, LIA_idx),[],2);

Ovx.(grupo).(animal).Spec.S_theta.or.SD= nanstd(S.or(:, theta_idx),[],2);
Ovx.(grupo).(animal).Spec.S_theta.pyr.SD= nanstd(S.pyr(:, theta_idx),[],2);
Ovx.(grupo).(animal).Spec.S_theta.rad.SD= nanstd(S.rad(:, theta_idx),[],2);
Ovx.(grupo).(animal).Spec.S_theta.slm.SD= nanstd(S.slm(:, theta_idx),[],2);
Ovx.(grupo).(animal).Spec.S_theta.hil.SD= nanstd(S.hil(:, theta_idx),[],2);


Ovx.(grupo).(animal).Spec.f= f;
Ovx.(grupo).(animal).Spec.t= t;
Ovx.(grupo).(animal).duracion.original= size(data,1);
Ovx.(grupo).(animal).duracion.ds= size(data_ds,1);
Ovx.(grupo).(animal).time_LIA_theta=[t_LIA, t_theta];

clearvars -except Ovx

cd 'E:\Users\database'
save ('ovario.mat', 'Ovx', '-append');

% save ('ovario.mat', 'Ovx', '-v7.3');


%% ----------------- 12) Plots comparativos
% genera plot comparativos
figLFP(Ovx);

figLFPnew_final(Ovx.Control, Ovx.Tratado);
%% ---------------- 13) Plot ejemplos
% Carga el archivo .abf deseado
list_abf = ls('*.abf*');             n=1;   % n= numero de abf en la lista         
data = abfload (strcat(pwd, '\', list_abf(n,:)));           
Fs= 20000;
t=(0:1/Fs:(size(data,1)-1)/Fs)';

% Cambiar intervalos segun la region que quieras plotear
int1 = 209*Fs; int2= 211.5*Fs; 
int1 = 124*Fs; int2= 126.5*Fs;
c= [0, 0, 0]; % Color en RGB

figure
hold off
%subplot(3,1,1)  % Pyramidal
plot(t(int1: int2), data(int1: int2, ch.pyr), 'Color', c); hold on
%subplot(3,1,2)  % Radiado
plot(t(int1: int2), data(int1: int2, ch.rad)-0.4, 'Color', c); hold on
%subplot(3,1,3)  % Slm
plot(t(int1: int2), data(int1: int2, ch.slm)-0.8, 'Color', c); hold on
%subplot(3,1,3)  % Giro dentado
%plot(t(int1: int2), data(int1: int2, ch.hil)-1.2, 'Color', c); hold on

% Para pintar ejes
plot([t(int1+0.2*Fs), t(int1+0.7*Fs)], [-1.2, -1.2], 'k')  % Eje X
text(t(int1+0.45*Fs), -1.25, '0.5 s')
plot([t(int1+0.2*Fs), t(int1+0.2*Fs)], [-1.2, +0.1-1.2], 'k')  % Eje Y
text(t(int1+0.1*Fs), -1.15, '0.1 mV')
xlim([t(int1), t(int2)]); ylim([-1.3, 0.4])
set(gca,'Xtick', [], 'Ytick', [], 'Box', 'off', 'FontSize', 14)