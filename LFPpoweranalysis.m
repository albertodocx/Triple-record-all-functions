
function [] = LFPpoweranalysis()

% Select *.mat file:

uiopen()

% Load variables:

LFPdata = results.LFP.LFPdata;
fs = results.LFP.fs;

% LFP profile:

pyrchannel = LFPdata(4, :);
fsdown = fs./1000;
newfs = 1000;
pyrch = downsample(pyrchannel, fsdown);

slmchannel = LFPdata(16, :);
slmch = downsample(slmchannel, fsdown);

pyrdata.dataLFP = pyrch;
pyrdata.fsLFP = newfs;

slmdata.dataLFP = slmch;
slmdata.fsLFP = newfs;

tWin = 2;
dtWin = 0.5;

pyrprofile = LFP_compute_psProfile(pyrdata, 'tWin', tWin, 'dtWin', dtWin);
slmprofile = LFP_compute_psProfile(slmdata, 'tWin', tWin, 'dtWin', dtWin);

%% PRE-PROCESAMIENTO DE LA SEÑAL (FILTROS, QUITAR RUIDO, ETC)



%% EXTRACT MEASURES + PLOTS

% Extract freq vector and power data for each channel:

freqspyr = pyrprofile.pSfreqs;
pspyr = pyrprofile.pS;

freqsslm = slmprofile.pSfreqs;
psslm = slmprofile.pS;

% Índices de fq de interés (valen para ambos canales)

thetaidx = find(freqspyr > 4 & freqspyr <= 12);
gammalowidx = find(freqspyr > 30 & freqspyr <= 60);
gammahighidx = find(freqspyr > 60 & freqspyr <= 90);
rippleidx = find(freqspyr > 150 & freqspyr <= 200);

% Extraigo power de ese rango para cada canal:
% Piramidal:

pyrchresults.thetapower = pspyr(thetaidx, :);
pyrchresults.gammalowpower = pspyr(gammalowidx, :);
pyrchresults.gammahighpower = pspyr(gammahighidx, :);
pyrchresults.ripplepower = pspyr(rippleidx, :);

% SLM:
slmchresults.thetapower = psslm(thetaidx, :);
slmchresults.gammalowpower = psslm(gammalowidx, :);
slmchresults.gammahighpower = psslm(gammahighidx, :);
slmchresults.ripplepower = psslm(rippleidx, :);


% Generar var de tiempo para power:

specfs = 1/dtWin;
time_spec = (0:(size(pspyr,2)-1))/specfs;

% Sacar mean por col de cada punto del power:

pyrchresults.meantheta = mean(pyrchresults.thetapower, 1); pyrchresults.meangammalow = mean(pyrchresults.gammalowpower, 1); pyrchresults.meangammahigh = mean(pyrchresults.gammahighpower, 1); pyrchresults.meanripple = mean(pyrchresults.ripplepower, 1);
slmchresults.meantheta = mean(slmchresults.thetapower, 1); slmchresults.meangammalow = mean(slmchresults.gammalowpower, 1); slmchresults.meangammahigh = mean(slmchresults.gammahighpower, 1); slmchresults.meanripple = mean(slmchresults.ripplepower, 1);

pyrchresults.meanthetasmoothed = smooth(pyrchresults.meantheta, 10, 'moving');pyrchresults.meangammalowsmoothed = smooth(pyrchresults.meangammalow, 10, 'moving');pyrchresults.meangammahighsmoothed = smooth(pyrchresults.meangammahigh, 10, 'moving'); pyrchresults.meanripplesmoothed = smooth(pyrchresults.meanripple, 10 , 'moving');
slmchresults.meanthetasmoothed = smooth(slmchresults.meantheta, 10, 'moving');slmchresults.meangammalowsmoothed = smooth(slmchresults.meangammalow, 10, 'moving');slmchresults.meangammahighsmoothed = smooth(slmchresults.meangammahigh, 10, 'moving'); slmchresults.meanripplesmoothed = smooth(slmchresults.meanripple, 10, 'moving');


% Espectrogramas:

climtheta = [0 1e5];
climgamlow = [0 1e4];
climgamhigh = [0 1e3];
climripple = [0 1e2] ;

fig3 = figure(3)
tiledlayout(3,1)

nexttile
%subplot(6, 1, 1)
imagesc(time_spec, freqspyr(thetaidx), pyrchresults.thetapower, climtheta)
title('Mean Power of Theta (4-10 Hz)')
colorbar
set(gca, 'YDir', 'normal')

nexttile
%subplot(6, 1, 2)
imagesc(time_spec, freqspyr(gammalowidx), pyrchresults.gammalowpower, climgamlow)
title('Mean Power of Gamma Low (30-60 Hz)')
colorbar
set(gca, 'YDir', 'normal')

nexttile
%subplot(6, 1, 3)
imagesc(time_spec, freqspyr(gammahighidx), pyrchresults.gammahighpower, climgamhigh)
title('Mean Power of Gamma High (60-90 Hz)')
colorbar
set(gca, 'YDir', 'normal')

% Ripple in pyr vs ripple in slm:
fig8 = figure(8)
tiledlayout(2, 1)
nexttile
imagesc(time_spec, freqspyr(rippleidx), pyrchresults.ripplepower)
title('Ripple Power (150- 200 Hz)')
colorbar
set(gca, 'YDir', 'normal')
nexttile
imagesc(time_spec, freqsslm(rippleidx), slmchresults.ripplepower)
title('Ripple Power SLM (150-200 Hz)')
colorbar
set(gca, 'YDir', 'normal')

% Theta vs Ripple in pyr:
fig9 = figure(9)
tiledlayout(2, 1)
nexttile
imagesc(time_spec, freqspyr(thetaidx), pyrchresults.thetapower, climtheta)
title('Mean Power of Theta (4-10 Hz)')
colorbar
set(gca, 'YDir', 'normal')
nexttile
imagesc(time_spec, freqspyr(rippleidx), pyrchresults.ripplepower, climripple)
title('Ripple Power (150- 200 Hz)')
colorbar
set(gca, 'YDir', 'normal')


fig4 = figure(4)
tiledlayout(3,1)

nexttile
%subplot(6, 1, 1)
imagesc(time_spec, freqsslm(thetaidx), slmchresults.thetapower, climtheta)
title('Mean Power of Theta (4-10 Hz)')
colorbar
set(gca, 'YDir', 'normal')

nexttile
%subplot(6, 1, 2)
imagesc(time_spec, freqsslm(gammalowidx), slmchresults.gammalowpower, climgamlow)
title('Mean Power of Gamma Low (30-60 Hz)')
colorbar
set(gca, 'YDir', 'normal')

nexttile
%subplot(6, 1, 3)
imagesc(time_spec, freqsslm(gammahighidx), slmchresults.gammahighpower, climgamhigh)
title('Mean Power of Gamma High (60-90 Hz)')
colorbar
set(gca, 'YDir', 'normal')

%% Selección manual de eventos de movimiento en el registro:

prompt = {'How many movement events do you detect?:'};

dlgtitle = 'Sample Buffer';
dims = [1 35];
definput = {'5'};

SB = inputdlg(prompt, dlgtitle, dims, definput); 

eventnr = str2num(SB{1});

for i = 1:eventnr

    %uiwait(msgbox('Select first left side, then right side of event', 'Instructions', "modal")); 
    [x, y, button] = ginput(2);

    leftidx = floor((x(1)).*1/dtWin); % dtWin porque es la fs para sacar el power (0,5 seg de desplazamiento por cada punto)  
    rightidx = round((x(2)).*1/dtWin); 

    %movevent = pyrchresults.thetapower(:, leftidx:rightidx); % Así tengo los valores de ese intervalo
    movevents.timestamps(i, :) = [time_spec(leftidx), time_spec(rightidx)];
    movevents.idx(i, :) = [leftidx, rightidx];

end

%% POWER COMO SEÑAL 

% Power como señal para comparar:

fig5 = figure(5)

tiledlayout(3, 1)

nexttile
plot(time_spec, pyrchresults.meantheta)
title('Mean Theta Power Pyramidal Channel')
nexttile
plot(time_spec, pyrchresults.meangammalow)
title('Mean Gamma Low Power Pyramidal Channel')
nexttile
plot(time_spec, pyrchresults.meangammahigh)
title('Mean Gamma High Power Pyramidal Channel')

fig6 = figure(6)

tiledlayout(3,1)

nexttile
plot(time_spec, slmchresults.meantheta)
title('Mean Theta Power SLM Channel')
nexttile
plot(time_spec, slmchresults.meangammalow)
title('Mean Gamma Low Power SLM Channel')
nexttile
plot(time_spec, slmchresults.meangammahigh)
title('Mean Gamma High Power SLM Channel')




end