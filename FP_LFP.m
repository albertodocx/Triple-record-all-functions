
function [] = FP_LFP()

% Function:

% Análisis LFP + FF + pos rueda. 

% 1. Load data:

% 1.1. Select TDT data:

uiwait(msgbox('Select "TDTdata.m" file', 'Instructions', "modal")); 

uiopen()

pyr = TDTdata.streams.Fie1.data.';
pyrfs = TDTdata.streams.Fie1.fs; % sólo haría falta una fs
gd = TDTdata.streams.Fie2.data.';
gdfs = TDTdata.streams.Fie2.fs;

% 1.2. Select results.m data:

uiwait(msgbox('Select results.m file', 'Instructions', "modal")); 

uiopen()

% Load variables: 

% 1.2.1. Cut LFP signal the same length as FP: 

SB = results.FP.params.SampleBuffer;
preSB = str2num(SB{1});
postSB = str2num(SB{2});

prepyr = round(preSB * pyrfs);
postpyr = round(postSB * pyrfs); 
pyr = pyr(prepyr+1:end-postpyr);

pregd = round(preSB * gdfs);
postgd = round(postSB * gdfs);
gd = gd(pregd +1: end-postgd);

% 1.2.2. Load bintime of wheel position (rows: wheelturns; cols = position
% from 1 to 10):

bintime = results.Analysis.wheelPosition.wheelnormcomplete.bintime;
wheelturns = size(bintime, 2);
posbin = length(fieldnames(bintime));

speed_threshold = results.Behavior.params.MovThreshold; 
wheel = results.Behavior.Position.processed;
speed = results.Behavior.Speed;

% 1.2.3. Get thresholded LFP data vs position:

diameter = pi * 40;
[pks, locs] = findpeaks(wheel, 'MinPeakHeight', diameter - 2);

posbin = [0:(diameter*0.1):diameter];

locs = [1, locs];
locs = [locs, length(pyr)]; 

pyr_thresh = pyr(speed >= speed_threshold);
gd_thresh = gd(speed >= speed_threshold);

idx_thresh = find(speed >= speed_threshold);
wheelval_thresh = wheel(speed >= speed_threshold);

newsigvectorwheel = zeros(length(wheel), 1).';
newsigvectorwheel(idx_thresh) = wheelval_thresh;
newsigvectorwheel(newsigvectorwheel == 0) = nan;


newsigvectorgd = zeros(length(gd), 1).';
newsigvectorgd(idx_thresh) = gd_thresh;
newsigvectorgd(newsigvectorgd == 0) = nan;

% 2. Plot LFP with speed threshold:

% Pyr channel:

fig1 = figure(1)
tiledlayout(3, 1)

nexttile
time = [0:length(pyr)-1]./pyrfs;
plot(time, speed, 'LineWidth', 0.001)
xlim([min(time) max(time)])
hold on
yline(speed_threshold, 'Color', 'red', 'LineWidth', 2)
ylabel('Speed (cm/s)')
hold off

nexttile

newsigvectorpyr = zeros(length(pyr), 1).';
newsigvectorpyr(idx_thresh) = pyr_thresh;
newsigvectorpyr(newsigvectorpyr == 0) = nan;

plot(time, newsigvectorpyr, 'LineWidth', 0.001)
xlim([min(time) max(time)])
hold on 
ylabel('Pyramidal Signal')
hold off


nexttile

newposvectorpyr = zeros(length(pyr), 1).';
newposvectorpyr(idx_thresh) = wheel(idx_thresh);
newposvectorpyr(newposvectorpyr == 0) = nan;

plot(time, newposvectorpyr, 'LineWidth', 0.001)
xlim([min(time) max(time)])
hold on
ylabel('Position (cm)')

xlabel('Time (s)')
hold off

% GD channel:

fig2 = figure(2)
tiledlayout(3, 1)

nexttile
time = [0:length(gd)-1]./gdfs;
plot(time, speed, 'LineWidth', 0.001)
xlim([min(time) max(time)])
hold on
yline(speed_threshold, 'Color', 'red', 'LineWidth', 2)
ylabel('Speed (cm/s)')
hold off

nexttile

newsigvectorgd = zeros(length(gd), 1).';
newsigvectorgd(idx_thresh) = gd_thresh;
newsigvectorgd(newsigvectorgd == 0) = nan;

plot(time, newsigvectorgd, 'LineWidth', 0.001)
xlim([min(time) max(time)])
hold on 
ylabel('Dentate Gyrus Signal')
hold off


nexttile

newposvectorgd = zeros(length(gd), 1).';
newposvectorgd(idx_thresh) = wheel(idx_thresh);
newposvectorgd(newposvectorgd == 0) = nan;

plot(time, newposvectorgd, 'LineWidth', 0.001)
xlim([min(time) max(time)])
hold on
ylabel('Position (cm)')

xlabel('Time (s)')
hold off





% 2. Plot LFP signal against position: 

for n = 1:(length(locs)-1)
        % Recortar por cada vuelta de la rueda y sacar vectores:
        temp = newsigvectorwheel(locs(n):locs(n+1)-1); 
        temppyr = newsigvectorpyr(locs(n):locs(n+1)-1);
        tempgd = newsigvectorpyr(locs(n):locs(n+1)-1);

    for p = 1:(length(posbin)-1)

        posname = sprintf('Pos%i', p);
        sigwthreshpyr(n).(posname) = temppyr((temp >= posbin(p)) & (temp < posbin(p+1)));
        sigwthreshgd(n).(posname) = tempgd((temp >= posbin(p)) & (temp < posbin(p+1)));
    
        bintimepyr(n).(posname) = [p-1:1./(length(sigwthreshpyr(n).(posname))-1):p];
        bintimegd(n).(posname) = [p-1:1./(length(sigwthreshgd(n).(posname))-1):p];

        % Calculate AUC:

        if length(bintimepyr(n).(posname)) > 200
            AUCdatapyr(n, p) = trapz(bintimepyr(n).(posname), sigwthreshpyr(n).(posname));    
            denspyr(n, p) = AUCdatapyr(n, p)./(length(sigwthreshpyr(n).(posname))./pyrfs);
            AUCdatagd(n, p) = trapz(bintimegd(n).(posname), sigwthreshgd(n).(posname));    
            densgd(n, p) = AUCdatagd(n, p)./(length(sigwthreshgd(n).(posname))./gdfs);
        else 
            AUCdatapyr(n, p) = nan;   
            denspyr(n, p) = nan;
            AUCdatagd(n, p) = nan;   
            densgd(n, p) = nan;
        end
        
    end

end

fig3 = figure(3)
ylim([-0.2 0.2])
hold on
for n = 1:(length(locs)-1)
    for p = 1:(length(posbin)-1)
        posname = sprintf('Pos%i', p);
        plot(bintimepyr(n).(posname), sigwthreshpyr(n).(posname), 'Color', 'black')
        hold on       
    end
end
ylabel('Pyramidal Channel')
hold off


fig4 = figure(4)
ylim([-0.2 0.2])
hold on
for n = 1:(length(locs)-1)
    for p = 1:(length(posbin)-1)
        posname = sprintf('Pos%i', p);
        plot(bintimegd(n).(posname), sigwthreshgd(n).(posname), 'Color', 'black')
        hold on       
    end
end
ylabel('DG Channel')
hold off

% Sin el threshold de velocidad: 



for n = 1:(length(locs)-1)
        % Recortar por cada vuelta de la rueda y sacar vectores:
        temp = wheel(locs(n):locs(n+1)-1); 
        temppyr = pyr(locs(n):locs(n+1)-1);
        tempgd = gd(locs(n):locs(n+1)-1);

    for p = 1:(length(posbin)-1)

        posname = sprintf('Pos%i', p);
        nothreshpyr(n).(posname) = temppyr((temp >= posbin(p)) & (temp < posbin(p+1)));
        nothreshgd(n).(posname) = tempgd((temp >= posbin(p)) & (temp < posbin(p+1)));
    
        bintimepyr(n).(posname) = [p-1:1./(length(nothreshpyr(n).(posname))-1):p];
        bintimegd(n).(posname) = [p-1:1./(length(nothreshgd(n).(posname))-1):p];

        % Calculate AUC:

        if length(bintimepyr(n).(posname)) > 200
            AUCdatapyr(n, p) = trapz(bintimepyr(n).(posname), nothreshpyr(n).(posname));    
            denspyr(n, p) = AUCdatapyr(n, p)./(length(nothreshpyr(n).(posname))./pyrfs);
            AUCdatagd(n, p) = trapz(bintimegd(n).(posname), nothreshgd(n).(posname));    
            densgd(n, p) = AUCdatagd(n, p)./(length(nothreshgd(n).(posname))./gdfs);
        else 
            AUCdatapyr(n, p) = nan;   
            denspyr(n, p) = nan;
            AUCdatagd(n, p) = nan;   
            densgd(n, p) = nan;
        end
        
    end

end


% Plots:


fig5 = figure(5)
ylim([-0.2 0.2])
hold on
for n = 1:(length(locs)-1)
    for p = 1:(length(posbin)-1)
        posname = sprintf('Pos%i', p);
        plot(bintimepyr(n).(posname), nothreshpyr(n).(posname), 'Color', 'black')
        hold on       
    end
end
ylabel('Pyramidal Channel')
hold off


fig6 = figure(6)
ylim([-0.5 0.5])
hold on
for n = 1:(length(locs)-1)
    for p = 1:(length(posbin)-1)
        posname = sprintf('Pos%i', p);
        plot(bintimegd(n).(posname), nothreshgd(n).(posname), 'Color', 'black')
        hold on       
    end
end
ylabel('DG Channel')
hold off


%% Análisis de la potencia/espectrograma:

% Parámetros (son los default ya en la fx)

tWin = 2;
dtWin = 0.5;

% Análisis potencia:

pyrdata.dataLFP = pyr;
pyrdata.fsLFP = pyrfs;

pyrch_power = LFP_compute_psProfile(pyrdata, 'tWin', tWin, 'dtWin', dtWin);

gddata.dataLFP = gd;
gddata.fsLFP = gdfs;

gdch_power = LFP_compute_psProfile(gddata, 'tWin', tWin, 'dtWin', dtWin);


%% Power de cada oscilación con rueda y FF:

% 1. Extraer vars:
% 
% 1.1. Señal FF:

GCAMP = results.FP.Signals.raw.GCAMP;
DFFZ = results.FP.Signals.DFFModZscore;

% 1.2. Power (como una señal): 

% De piramidal y giro dentado (de los dos canales que hay guardados):

freqspyr = pyrch_power.pSfreqs;
pspyr = pyrch_power.pS;

freqsgd = gdch_power.pSfreqs;
psgd = gdch_power.pS;

% Índices de fq de interés (valen para ambos canales)

thetaidx = find(freqspyr > 4 & freqspyr <= 10);
gammalowidx = find(freqspyr > 30 & freqspyr <= 60);
gammahighidx = find(freqspyr > 60 & freqspyr <= 90);

% Extraigo power de ese rango para cada canal:
% Piramidal:

pyrch.thetapower = pspyr(thetaidx, :);
pyrch.gammalowpower = pspyr(gammalowidx, :);
pyrch.gammahighpower = pspyr(gammahighidx, :);

% Giro dentado:
gdch.thetapower = psgd(thetaidx, :);
gdch.gammalowpower = psgd(gammalowidx, :);
gdch.gammahighpower = psgd(gammahighidx, :);

% 2. Generar plots de power de cada oscilación con la rueda y la FF:

% Generar var de tiempo para power y para ff:

specfs = 1/dtWin;
time_spec = (0:(size(pspyr,2)-1))/specfs;

fsFP = results.FP.params.fs;
time_FP = [0:length(GCAMP)-1]./fsFP;

fswheel = results.Behavior.Fs;
time_wheel = [0:length(wheel)-1]./fswheel;

% Sacar mean por col de cada punto del power:

pyrch.meantheta = mean(pyrch.thetapower, 1); pyrch.meangammalow = mean(pyrch.gammalowpower, 1); pyrch.meangammahigh = mean(pyrch.gammahighpower, 1);
gdch.meantheta = mean(gdch.thetapower, 1); gdch.meangammalow = mean(gdch.gammalowpower, 1); gdch.meangammahigh = mean(gdch.gammahighpower, 1);

pyrch.meanthetasmoothed = smooth(pyrch.meantheta, 10, 'moving');pyrch.meangammalowsmoothed = smooth(pyrch.meangammalow, 10, 'moving');pyrch.meangammahighsmoothed = smooth(pyrch.meangammahigh, 10, 'moving');
gdch.meanthetasmoothed = smooth(gdch.meantheta, 10, 'moving');gdch.meangammalowsmoothed = smooth(gdch.meangammalow, 10, 'moving');gdch.meangammahighsmoothed = smooth(gdch.meangammahigh, 10, 'moving');



% Plot de FF + potencia + rueda de canal piramidal:

fig7 = figure('Name', 'Pyramidal Channel')
subplot(6, 1, 1)
plot(time_FP, GCAMP)
title('Raw GCAMP')

subplot(6, 1, 2)
plot(time_FP, DFFZ)
title('Df/f Z-Score FP')

subplot(6, 1, 3)
plot(time_wheel, wheel)
title('Wheel Signal')

subplot(6, 1, 4)
plot(time_spec, pyrch.meantheta)
title('Mean Power of Theta (4-10 Hz)')

subplot(6, 1, 5)
plot(time_spec, pyrch.meangammalow)
title('Mean Power of Gamma Low (30-60 Hz)')

subplot(6, 1, 6)
plot(time_spec, pyrch.meangammahigh)
title('Mean Power of Gamma High (60-90 Hz)')


% De canal de GD:

fig8 = figure('Name', 'Dentate Gyrus Channel')
subplot(6, 1, 1)
plot(time_FP, GCAMP)
title('Raw GCAMP')

subplot(6, 1, 2)
plot(time_FP, DFFZ)
title('Df/f Z-Score FP')

subplot(6, 1, 3)
plot(time_wheel, wheel)
title('Wheel Signal')

subplot(6, 1, 4)
plot(time_spec, gdch.meantheta)
title('Mean Power of Theta (4-10 Hz)')

subplot(6, 1, 5)
plot(time_spec, gdch.meangammalow)
title('Mean Power of Gamma Low (30-60 Hz)')

subplot(6, 1, 6)
plot(time_spec, gdch.meangammahigh)
title('Mean Power of Gamma High (60-90 Hz)')

% Con Smooth
% Pyr:
fig9 = figure('Name', 'Pyramidal Channel - LFP smoothed')
subplot(6, 1, 1)
plot(time_FP, GCAMP)
title('Raw GCAMP')

subplot(6, 1, 2)
plot(time_FP, DFFZ)
title('Df/f Z-Score FP')

subplot(6, 1, 3)
plot(time_wheel, wheel)
title('Wheel Signal')

subplot(6, 1, 4)
plot(time_spec, pyrch.meanthetasmoothed)
title('Mean Power of Theta (4-10 Hz)')

subplot(6, 1, 5)
plot(time_spec, pyrch.meangammalowsmoothed)
title('Mean Power of Gamma Low (30-60 Hz)')

subplot(6, 1, 6)
plot(time_spec, pyrch.meangammahighsmoothed)
title('Mean Power of Gamma High (60-90 Hz)')

%GD:
fig10 = figure('Name', 'Dentate Gyrus Channel - LFP smoothed')
subplot(6, 1, 1)
plot(time_FP, GCAMP)
title('Raw GCAMP')

subplot(6, 1, 2)
plot(time_FP, DFFZ)
title('Df/f Z-Score FP')

subplot(6, 1, 3)
plot(time_wheel, wheel)
title('Wheel Signal')

subplot(6, 1, 4)
plot(time_spec, gdch.meanthetasmoothed)
title('Mean Power of Theta (4-10 Hz)')

subplot(6, 1, 5)
plot(time_spec, gdch.meangammalowsmoothed)
title('Mean Power of Gamma Low (30-60 Hz)')

subplot(6, 1, 6)
plot(time_spec, gdch.meangammahighsmoothed)
title('Mean Power of Gamma High (60-90 Hz)')


% 3. Espectogramas con rangos de frecuencias en el tiempo:

% Pyr
clims = [0 1e-4];
fig11 = figure('Name', 'Spectrogram Pyramidal Channel')
subplot(6, 1, 1)
plot(time_FP, GCAMP)
title('Raw GCAMP')

subplot(6, 1, 2)
plot(time_FP, DFFZ)
title('Df/f Z-Score FP')

subplot(6, 1, 3)
plot(time_wheel, wheel)
title('Wheel Signal')

subplot(6, 1, 4)
imagesc(time_spec, freqspyr(thetaidx), pyrch.thetapower, clims)
title('Mean Power of Theta (4-10 Hz)')
set(gca, 'YDir', 'normal')

subplot(6, 1, 5)
imagesc(time_spec, freqspyr(gammalowidx), pyrch.gammalowpower, clims)
title('Mean Power of Gamma Low (30-60 Hz)')
set(gca, 'YDir', 'normal')

subplot(6, 1, 6)
imagesc(time_spec, freqspyr(gammahighidx), pyrch.gammahighpower, clims)
title('Mean Power of Gamma High (60-90 Hz)')
set(gca, 'YDir', 'normal')


% GD

fig12 = figure('Name', 'Spectrogram Dentate Gyrus Channel')
subplot(6, 1, 1)
plot(time_FP, GCAMP)
title('Raw GCAMP')

subplot(6, 1, 2)
plot(time_FP, DFFZ)
title('Df/f Z-Score FP')

subplot(6, 1, 3)
plot(time_wheel, wheel)
title('Wheel Signal')

subplot(6, 1, 4)
imagesc(time_spec, freqspyr(thetaidx), gdch.thetapower, clims)
title('Mean Power of Theta (4-10 Hz)')
set(gca, 'YDir', 'normal')

subplot(6, 1, 5)
imagesc(time_spec, freqspyr(gammalowidx), gdch.gammalowpower, clims)
title('Mean Power of Gamma Low (30-60 Hz)')
set(gca, 'YDir', 'normal')

subplot(6, 1, 6)
imagesc(time_spec, freqspyr(gammahighidx), gdch.gammahighpower, clims)
title('Mean Power of Gamma High (60-90 Hz)')
set(gca, 'YDir', 'normal')

% 4. Oscilación en función de la posición

% CANAL PIRAMIDAL:

% Con locs en base a las muestras de la rueda puedo encontrar segundos en
% power y sacar la posición a partir de ahí:

locstime = zeros([1 length(locs)]);
locstime(2:end) = locs(2:end)/fswheel;

for n = 1:(length(locs)-1)
        % Recortar por cada vuelta de la rueda y sacar vectores:
        temp = wheel(locs(n):locs(n+1)-1); 
        temptimeidx = [find(time_spec >= locstime(n), 1, 'first'):find(time_spec < locstime(n+1),1, 'last')]; % La vuelta completa desde el pico hasta el siguiente pico o el final
        pyrtemptheta = pyrch.thetapower(:, temptimeidx);
        pyrtempgammalow = pyrch.gammalowpower(:, temptimeidx);
        pyrtempgammahigh = pyrch.gammahighpower(:, temptimeidx);
        gdtemptheta = gdch.thetapower(:, temptimeidx);
        gdtempgammalow = gdch.gammalowpower(:, temptimeidx);
        gdtempgammahigh = gdch.gammahighpower(:, temptimeidx);
        
        len_temp = length(temp);
        len_pow = length(temptimeidx);

    for p = 1:(length(posbin)-1)

        posname = sprintf('Pos%i', p);
        wheelidx = find( temp >= posbin(p) & temp < posbin(p+1)); % no está bien calculado porque te recorta el vector y ya no sabes qué índices tienes
        wheel_t = wheelidx/fswheel;
        t_idx = [1:length(pyrtemptheta)]/2; % Muestreo del power es de 2 Hz (2 muestras cada segundo o lo que es lo mismo una dtWin de 0,5)
        %t_idx = find(t_idx >= wheel_t(1), 1, 'first'):find(t_idx <= wheel_t(end), 1, 'last');
        if ~isempty(wheel_t)

        t_idx = [find(t_idx >= wheel_t(1), 1, 'first'):find(t_idx < wheel_t(end), 1, 'last')];

% Canal piramidal
        pyrpos.postheta(n).(posname) = pyrtemptheta(:, t_idx);
        pyrpos.posgammalow(n).(posname) = pyrtempgammalow(:, t_idx);
        pyrpos.posgammahigh(n).(posname) = pyrtempgammahigh(:, t_idx);
        pyrpos.meanpostheta(n).(posname) = mean(pyrpos.postheta(n).(posname), 1);
        pyrpos.meanposgammalow(n).(posname) = mean(pyrpos.posgammalow(n).(posname), 1);
        pyrpos.meanposgammahigh(n).(posname) = mean(pyrpos.posgammahigh(n).(posname), 1);
    
        pyrpos.bintimetheta(n).(posname) = [p-1:1./(size(pyrpos.postheta(n).(posname), 2)-1):p];
        pyrpos.bintimegammalow(n).(posname) = [p-1:1./(size(pyrpos.posgammalow(n).(posname), 2)-1):p];
        pyrpos.bintimegammahigh(n).(posname) = [p-1:1./(size(pyrpos.posgammahigh(n).(posname), 2)-1):p];
% Canal giro dentado
        gdpos.postheta(n).(posname) = gdtemptheta(:, t_idx);
        gdpos.posgammalow(n).(posname) = gdtempgammalow(:, t_idx);
        gdpos.posgammahigh(n).(posname) = gdtempgammahigh(:, t_idx);
        gdpos.meanpostheta(n).(posname) = mean(gdpos.postheta(n).(posname), 1);
        gdpos.meanposgammalow(n).(posname) = mean(gdpos.posgammalow(n).(posname), 1);
        gdpos.meanposgammahigh(n).(posname) = mean(gdpos.posgammahigh(n).(posname), 1);
    
        gdpos.bintimetheta(n).(posname) = [p-1:1./(size(gdpos.postheta(n).(posname), 2)-1):p];
        gdpos.bintimegammalow(n).(posname) = [p-1:1./(size(gdpos.posgammalow(n).(posname), 2)-1):p];
        gdpos.bintimegammahigh(n).(posname) = [p-1:1./(size(gdpos.posgammahigh(n).(posname), 2)-1):p];

        
        else
        
        pyrpos.postheta(n).(posname) = [];
        pyrpos.posgammalow(n).(posname) = [];
        pyrpos.posgammahigh(n).(posname) = [];
        pyrpos.meanpostheta(n).(posname) = [];
        pyrpos.meanposgammalow(n).(posname) = [];
        pyrpos.meanposgammahigh(n).(posname) = [];
    
        pyrpos.bintimetheta(n).(posname) = [];
        pyrpos.bintimegammalow(n).(posname) = [];
        pyrpos.bintimegammahigh(n).(posname) = [];


        gdpos.postheta(n).(posname) = [];
        gdpos.posgammalow(n).(posname) = [];
        gdpos.posgammahigh(n).(posname) = [];
        gdpos.meanpostheta(n).(posname) = [];
        gdpos.meanposgammalow(n).(posname) = [];
        gdpos.meanposgammahigh(n).(posname) = [];
    
        gdpos.bintimetheta(n).(posname) = [];
        gdpos.bintimegammalow(n).(posname) = [];
        gdpos.bintimegammahigh(n).(posname) = [];


        end

        if ~isempty(pyrpos.bintimetheta(n).(posname))

            pyrpos.AUCtheta(n, p) = trapz(pyrpos.bintimetheta(n).(posname), pyrpos.meanpostheta(n).(posname));
            pyrpos.AUCgammalow(n, p) = trapz(pyrpos.bintimegammalow(n).(posname), pyrpos.meanposgammalow(n).(posname));
            pyrpos.AUCgammahi(n, p) = trapz(pyrpos.bintimegammahigh(n).(posname), pyrpos.meanposgammahigh(n).(posname));
    
            pyrpos.denstheta(n, p) = pyrpos.AUCtheta(n, p)./(length(pyrpos.meanpostheta(n).(posname))./2); % 2 Hz
            pyrpos.densgammalow(n, p) = pyrpos.AUCgammalow(n, p)./(length(pyrpos.meanpostheta(n).(posname))./2);
            pyrpos.densgammahi(n, p) = pyrpos.AUCgammahi(n, p)./(length(pyrpos.meanpostheta(n).(posname))./2);

            gdpos.AUCtheta(n, p) = trapz(gdpos.bintimetheta(n).(posname), gdpos.meanpostheta(n).(posname));
            gdpos.AUCgammalow(n, p) = trapz(gdpos.bintimegammalow(n).(posname), gdpos.meanposgammalow(n).(posname));
            gdpos.AUCgammahi(n, p) = trapz(gdpos.bintimegammahigh(n).(posname), gdpos.meanposgammahigh(n).(posname));
    
            gdpos.denstheta(n, p) = gdpos.AUCtheta(n, p)./(length(gdpos.meanpostheta(n).(posname))./2); % 2 Hz
            gdpos.densgammalow(n, p) = gdpos.AUCgammalow(n, p)./(length(gdpos.meanposgammalow(n).(posname))./2);
            gdpos.densgammahi(n, p) = gdpos.AUCgammahi(n, p)./(length(gdpos.meanposgammahigh(n).(posname))./2);

        else 
            pyrpos.AUCtheta(n, p) = nan;
            pyrpos.AUCgammalow(n, p) = nan;
            pyrpos.AUCgammahi(n, p) = nan;
    
            pyrpos.denstheta(n, p) = nan; % 2 Hz
            pyrpos.densgammalow(n, p) = nan;
            pyrpos.densgammahi(n, p) = nan;

            gdpos.AUCtheta(n, p) = nan;
            gdpos.AUCgammalow(n, p) = nan;
            gdpos.AUCgammahi(n, p) = nan;
    
            gdpos.denstheta(n, p) = nan; % 2 Hz
            gdpos.densgammalow(n, p) = nan;
            gdpos.densgammahi(n, p) = nan;
        end
    end

end


% Plots en función de la posición:

% Piramidal:
fig13 = figure('Name', 'Pyr Channel. PowerXPosition')
subplot(3, 1, 1)
%ylim([0 1e-4])
hold on
for n = 1:(length(locs)-1)
    for p = 1:(length(posbin)-1)
        posname = sprintf('Pos%i', p);
        plot(pyrpos.bintimetheta(n).(posname), pyrpos.meanpostheta(n).(posname), 'Color', 'black')
        hold on       
    end
end
title('Mean Theta Power - Pyramidal Channel')
hold off


subplot(3, 1, 2)
%ylim([0 1e-5])
hold on
for n = 1:(length(locs)-1)
    for p = 1:(length(posbin)-1)
        posname = sprintf('Pos%i', p);
        plot(pyrpos.bintimegammalow(n).(posname), pyrpos.meanposgammalow(n).(posname), 'Color', 'black')
        hold on       
    end
end
title('Mean Gamma Low Power - Pyramidal Channel')
hold off

subplot(3, 1, 3)
%ylim([0 1e-6])
hold on
for n = 1:(length(locs)-1)
    for p = 1:(length(posbin)-1)
        posname = sprintf('Pos%i', p);
        plot(pyrpos.bintimegammahigh(n).(posname), pyrpos.meanposgammahigh(n).(posname), 'Color', 'black')
        hold on       
    end
end
title('Mean Gamma High Power - Pyramidal Channel')
hold off



% Giro dentado:

fig14 = figure('Name', 'GD Channel. PowerXPosition')
subplot(3, 1, 1)
%ylim([0 1e-4])
hold on
for n = 1:(length(locs)-1)
    for p = 1:(length(posbin)-1)
        posname = sprintf('Pos%i', p);
        plot(gdpos.bintimetheta(n).(posname), gdpos.meanpostheta(n).(posname), 'Color', 'black')
        hold on       
    end
end
title('Mean Theta Power - GD Channel')
hold off


subplot(3, 1, 2)
%ylim([0 1e-5])
hold on
for n = 1:(length(locs)-1)
    for p = 1:(length(posbin)-1)
        posname = sprintf('Pos%i', p);
        plot(gdpos.bintimegammalow(n).(posname), gdpos.meanposgammalow(n).(posname), 'Color', 'black')
        hold on       
    end
end
title('Mean Gamma Low Power - GD Channel')
hold off

subplot(3, 1, 3)
%ylim([0 1e-6])
hold on
for n = 1:(length(locs)-1)
    for p = 1:(length(posbin)-1)
        posname = sprintf('Pos%i', p);
        plot(gdpos.bintimegammahigh(n).(posname), gdpos.meanposgammahigh(n).(posname), 'Color', 'black')
        hold on       
    end
end
title('Mean Gamma High Power - GD Channel')
hold off


% 4. Espectograma donde cada columna es un bin y el eje y es un espectro de
% frecuencias: 


for p = 1:(length(posbin)-1)
    posname = sprintf('Pos%i', p);
    for n = 1:(length(locs)-1)
        if ~isempty(pyrpos.postheta(n).(posname))
            posmean_pyrtheta(:, n) = mean(pyrpos.postheta(n).(posname), 2);
            posmean_pyrgamlow(:, n) = mean(pyrpos.posgammalow(n).(posname), 2);
            posmean_pyrgamhi(:, n) = mean(pyrpos.posgammahigh(n).(posname), 2);

            posmean_gdtheta(:, n) = mean(gdpos.postheta(n).(posname), 2);
            posmean_gdgamlow(:, n) = mean(gdpos.posgammalow(n).(posname), 2);
            posmean_gdgamhi(:, n) = mean(gdpos.posgammahigh(n).(posname), 2);
        else
            posmean_pyrtheta(:, n) = nan([length(thetaidx) 1]);
            posmean_pyrgamlow(:, n) = nan([length(gammalowidx) 1]);
            posmean_pyrgamhi(:, n) = nan([length(gammahighidx) 1]);

            posmean_gdtheta(:, n) = nan([length(thetaidx) 1]);
            posmean_gdgamlow(:, n) = nan([length(gammalowidx) 1]);
            posmean_gdgamhi(:, n) = nan([length(gammahighidx) 1]);        
        end
    end
    turnmean_pyrtheta(:, p) = nanmean(posmean_pyrtheta, 2);
    turnmean_pyrgamlow(:, p) = nanmean(posmean_pyrgamlow, 2);
    turnmean_pyrgamhi(:, p) = nanmean(posmean_pyrgamhi, 2);

    turnmean_gdtheta(:, p) = nanmean(posmean_gdtheta, 2);
    turnmean_gdgamlow(:, p) = nanmean(posmean_gdgamlow, 2);
    turnmean_gdgamhi(:, p) = nanmean(posmean_gdgamhi, 2);
end


% Espectrogramas por posición:

xval = [0:1:10];
clim = [0 1e-5];

% Piramidal:
fig15 = figure('Name', 'Spectrogram: LFPXPosition Pyramidal')
subplot(3, 1, 1)
imagesc(xval, freqspyr(thetaidx), turnmean_pyrtheta, clim)
%contourf(turnmean_pyrtheta)
set(gca, 'YDir', 'normal')
title('Mean Theta Power - Pyr Channel')


subplot(3, 1, 2)
imagesc(xval, freqspyr(gammalowidx), turnmean_pyrgamlow, clim)
%contourf(turnmean_pyrgamlow)
set(gca, 'YDir', 'normal')
title('Mean Low Gamma Power - Pyr Channel')


subplot(3, 1, 3)
imagesc(xval, freqspyr(gammahighidx), turnmean_pyrgamhi, clim)
%contourf(turnmean_pyrgamhi)
set(gca, 'YDir', 'normal')
title('Mean High Gamma Power - Pyr Channel')

% GD: (muy feos,.. probar con contourf.. u otros similares)
fig16 = figure('Name', 'Spectrogram: LFPXPosition DG Channel')
subplot(3, 1, 1)
imagesc(xval, freqspyr(thetaidx), turnmean_gdtheta, clim)
%contourf(turnmean_pyrtheta)
set(gca, 'YDir', 'normal')
title('Mean Theta Power - GD Channel')


subplot(3, 1, 2)
imagesc(xval, freqspyr(gammalowidx), turnmean_gdgamlow, clim)
%contourf(turnmean_pyrgamlow)
set(gca, 'YDir', 'normal')
title('Mean Low Gamma Power - GD Channel')


subplot(3, 1, 3)
imagesc(xval, freqspyr(gammahighidx), turnmean_gdgamhi, clim)
%contourf(turnmean_pyrgamhi)
set(gca, 'YDir', 'normal')
title('Mean High Gamma Power - GD Channel')


% 5. ÁREA Y DENSIDAD DE POWER EN FUNCIÓN DE LA POSICIÓN

% AREA
pyrpos.AUCthetamean = nanmean(pyrpos.AUCtheta, 1).';
pyrpos.AUCgammalowmean = nanmean(pyrpos.AUCgammalow, 1).';
pyrpos.AUCgammahimean = nanmean(pyrpos.AUCgammahi, 1).';

posnames = fieldnames(pyrpos.bintimetheta);
pyrAUCdatatable = table(posnames, pyrpos.AUCthetamean, pyrpos.AUCgammalowmean, pyrpos.AUCgammahimean, 'VariableNames', ["Posnames" "Theta" "GammaLow" "GammaHigh"]);
pyrAUCdatatable.Posnames = categorical(pyrAUCdatatable.Posnames);


gdpos.AUCthetamean = nanmean(gdpos.AUCtheta, 1).';
gdpos.AUCgammalowmean = nanmean(gdpos.AUCgammalow, 1).';
gdpos.AUCgammahimean = nanmean(gdpos.AUCgammahi, 1).';

posnames = fieldnames(gdpos.bintimetheta);
gdAUCdatatable = table(posnames, gdpos.AUCthetamean, gdpos.AUCgammalowmean, gdpos.AUCgammahimean, 'VariableNames', ["Posnames" "Theta" "GammaLow" "GammaHigh"]);
gdAUCdatatable.Posnames = categorical(gdAUCdatatable.Posnames);

% DENSITY


pyrpos.densthetamean = nanmean(pyrpos.denstheta, 1).';
pyrpos.densgammalowmean = nanmean(pyrpos.densgammalow, 1).';
pyrpos.densgammahimean = nanmean(pyrpos.densgammahi, 1).';

posnames = fieldnames(pyrpos.bintimetheta);
pyrdensdatatable = table(posnames, pyrpos.densthetamean, pyrpos.densgammalowmean, pyrpos.densgammahimean, 'VariableNames', ["Posnames" "Theta" "GammaLow" "GammaHigh"]);
pyrdensdatatable.Posnames = categorical(pyrdensdatatable.Posnames);


gdpos.densthetamean = nanmean(gdpos.denstheta, 1).';
gdpos.densgammalowmean = nanmean(gdpos.densgammalow, 1).';
gdpos.densgammahimean = nanmean(gdpos.densgammahi, 1).';

posnames = fieldnames(gdpos.bintimetheta);
gddensdatatable = table(posnames, gdpos.densthetamean, gdpos.densgammalowmean, gdpos.densgammahimean, 'VariableNames', ["Posnames" "Theta" "GammaLow" "GammaHigh"]);
gddensdatatable.Posnames = categorical(gddensdatatable.Posnames);



% Bar-plots:

fig17 = figure(17)

tiledlayout(3, 1)

nexttile
bar(pyrAUCdatatable.Theta);
hold on 
ylabel('AUC mean - Theta')
hold off

nexttile
bar(pyrAUCdatatable.GammaLow);
hold on
ylabel('AUC mean - Gamma Low')
hold off


nexttile
bar(pyrAUCdatatable.GammaHigh);
hold on
ylabel('AUC mean - Gamma High')
hold off

fig18 = figure(18)
tiledlayout(3, 1)

nexttile
bar(gdAUCdatatable.Theta);
hold on 
ylabel('AUC mean - Theta')
hold off

nexttile
bar(gdAUCdatatable.GammaLow);
hold on
ylabel('AUC mean - Gamma Low')
hold off

nexttile
bar(gdAUCdatatable.GammaHigh);
hold on
ylabel('AUC mean - Gamma High')
hold off

% Density:

fig19 = figure(19)
tiledlayout(3, 1)

nexttile
bar(pyrdensdatatable.Theta);
hold on 
ylabel('Density mean - Theta')
hold off

nexttile
bar(pyrdensdatatable.GammaLow);
hold on
ylabel('Density mean - Gamma Low')
hold off

nexttile
bar(pyrdensdatatable.GammaHigh);
hold on
ylabel('Density mean - Gamma High')
hold off

fig20 = figure(20)

tiledlayout(3, 1)

nexttile
bar(gddensdatatable.Theta);
hold on 
ylabel('Density mean - Theta')
hold off

nexttile
bar(gddensdatatable.GammaLow);
hold on
ylabel('Density mean - Gamma Low')
hold off

nexttile
bar(gddensdatatable.GammaHigh);
hold on
ylabel('Density mean - Gamma High')
hold off


% Boxplots de power por cada pos:


tiledlayout(3, 1)

nexttile
boxchart(turnmean_pyrtheta)
hold on
ylabel("Theta Power")
title("Pyramidal Channel")
hold off

nexttile
boxchart(turnmean_pyrgamlow)
hold on
ylabel("Gamma Low Power")
hold off

nexttile
boxchart(turnmean_pyrgamhi)
hold on
ylabel("Gamma High Power")
hold off



tiledlayout(3, 1)

nexttile
boxchart(turnmean_gdtheta)
hold on
ylabel("Theta Power")
title("DG Channel")
hold off

nexttile
boxchart(turnmean_gdgamlow)
hold on
ylabel("Gamma Low Power")
hold off

nexttile
boxchart(turnmean_gdgamhi)
hold on
ylabel("Gamma High Power")
hold off





%%% Pendiente: lo mismo con GD y lo mismo para ambos canales con densidad

%% Find Theta

pyrdata.dataLFP = pyr;
pyrdata.fsLFP = pyrfs;

gddata.dataLFP = gd;
gddata.fsLFP = gdfs;

pyr_theta_analysis_peak = LFP_find_theta(pyrdata, 'peak_trough', 'peak');
pyr_theta_analysis_trough = LFP_find_theta(pyrdata, 'peak_trough', 'trough');

gd_theta_analysis_peak = LFP_find_theta(gddata, 'peak_trough', 'peak');
gd_theta_analysis_trough = LFP_find_theta(gddata, 'peak_trough', 'trough');

% 1. Sacar índices de picos y valles de theta:

idx_pyr_peak = pyr_theta_analysis_peak.iCyc;
idx_pyr_trough =  pyr_theta_analysis_trough.iCyc;

idx_gd_peak = gd_theta_analysis_peak.iCyc;
idx_gd_trough = gd_theta_analysis_trough.iCyc;

% 2. Plots con FF:

% Pyramidal
iCyc = idx_pyr_peak;
tiledlayout(2, 2)
nexttile
for i = 1:size(iCyc, 2)
    idx = iCyc(1, i):iCyc(2, i);
    time = [0:(length(idx)-1)]/pyrfs;
    pyrpeak_timestamps(i).time = time;
    values = smooth(pyr(idx), 100, 'moving');
    pyrpeak_signal(i).value = values;
    pyrpeak_longd(i) = length(pyrpeak_timestamps(i).time); % Guardar longitud
    plot(time, pyr(idx), 'black')
    hold on
end
ylabel('Peaks')
hold off


nexttile
for i = 1:size(iCyc, 2)
    idx = iCyc(1, i):iCyc(2, i);
    time = [0:(length(idx)-1)]/fsFP;
    pyrpeak_timeDFF(i).time = time;
    values = smooth(DFFZ(idx), 100, 'moving');
    pyrpeak_DFFsignal(i).value = values;
    plot(time, DFFZ(idx), 'black')
    hold on
end
ylabel('DFF Peaks')
hold off

iCyc = idx_pyr_trough;
nexttile
for i = 1:size(iCyc, 2)
    idx = iCyc(1, i):iCyc(2, i);
    time = [0:(length(idx)-1)]/pyrfs;
    pyrtrough_timestamps(i).time = time;
    values = smooth(pyr(idx), 100, 'moving');
    pyrtrough_signal(i).value = values;
    pyrtrough_longd(i) = length(pyrtrough_timestamps(i).time); 
    plot(time, pyr(idx), 'black')
    hold on
end
ylabel('Troughs')
hold off


nexttile
for i = 1:size(iCyc, 2)
    idx = iCyc(1, i):iCyc(2, i);
    time = [0:(length(idx)-1)]/fsFP;
    pyrtrough_timeDFF(i).time = time;
    values = smooth(DFFZ(idx), 100, 'moving');
    pyrtrough_DFFsignal(i).value = values;
    plot(time, DFFZ(idx), 'black')
    hold on
end
ylabel('DFF Troughs')
hold off

% giro dentado:

iCyc = idx_gd_peak;
tiledlayout(2, 2)
nexttile
for i = 1:size(iCyc, 2)
    idx = iCyc(1, i):iCyc(2, i);
    time = [0:(length(idx)-1)]/gdfs;
    gdpeak_timestamps(i).time = time;
    values = smooth(gd(idx), 100, 'moving');
    gdpeak_signal(i).value = values;
    gdpeak_longd(i) = length(gdpeak_timestamps(i).time); 
    plot(time, gd(idx), 'black')
    hold on
end
ylabel('Peaks')
hold off


nexttile
for i = 1:size(iCyc, 2)
    idx = iCyc(1, i):iCyc(2, i);
    time = [0:(length(idx)-1)]/fsFP;
    gdpeak_timeDFF(i).time = time;
    values = smooth(DFFZ(idx), 100, 'moving');
    gdpeak_DFFsignal(i).value = values;
    plot(time, DFFZ(idx), 'black')
    hold on
end
ylabel('DFF Peaks')
hold off

iCyc = idx_gd_trough;
nexttile
for i = 1:size(iCyc, 2)
    idx = iCyc(1, i):iCyc(2, i);
    time = [0:(length(idx)-1)]/gdfs;
    gdtrough_timestamps(i).time = time;
    values = smooth(gd(idx), 100, 'moving');
    gdtrough_signal(i).value = values;
    gdtrough_longd(i) = length(gdtrough_timestamps(i).time); 
    plot(time, gd(idx), 'black')
    hold on
end
ylabel('Troughs')
hold off


nexttile
for i = 1:size(iCyc, 2)
    idx = iCyc(1, i):iCyc(2, i);
    time = [0:(length(idx)-1)]/fsFP;
    gdtrough_timeDFF(i).time = time;
    values = smooth(DFFZ(idx), 100, 'moving');
    gdtrough_DFFsignal(i).value = values;
    plot(time, DFFZ(idx), 'black')
    hold on
end
ylabel('DFF Troughs')
hold off

% Reajustar las matrices para imagesc:

[val, ~] = max(pyrpeak_longd);

for i = 1:size(idx_pyr_peak, 2)
    if length(pyrpeak_DFFsignal(i).value) < val
       sam = val - length(pyrpeak_DFFsignal(i).value); 
       pyrpeak_timeDFF(i).time(end+1:end+sam) = length(pyrpeak_timeDFF(i).time) + [1:1:sam].';
       pyrpeak_DFFsignal(i).value(end+1:end+sam) = nan;
       pyrpeak_timestamps(i).time(end+1:end+sam) = length(pyrpeak_timestamps(i).time) + [1:1:sam].';
       pyrpeak_signal(i).value(end+1:end+sam) = nan;
    end
end

[val, ~] = max(pyrtrough_longd);

for i = 1:size(idx_pyr_trough, 2)
    if length(pyrtrough_DFFsignal(i).value) < val
       sam = val - length(pyrtrough_DFFsignal(i).value); 
       pyrtrough_timeDFF(i).time(end+1:end+sam) = length(pyrtrough_timeDFF(i).time) + [1:1:sam].';
       pyrtrough_DFFsignal(i).value(end+1:end+sam) = nan;
       pyrtrough_timestamps(i).time(end+1:end+sam) = length(pyrtrough_timestamps(i).time) + [1:1:sam].';
       pyrtrough_signal(i).value(end+1:end+sam) = nan;
    end
end


[val, ~] = max(gdpeak_longd);

for i = 1:size(idx_gd_peak, 2)
    if length(gdpeak_DFFsignal(i).value) < val
       sam = val - length(gdpeak_DFFsignal(i).value); 
       gdpeak_timeDFF(i).time(end+1:end+sam) = length(gdpeak_timeDFF(i).time) + [1:1:sam].';
       gdpeak_DFFsignal(i).value(end+1:end+sam) = nan;
       gdpeak_timestamps(i).time(end+1:end+sam) = length(gdpeak_timestamps(i).time) + [1:1:sam].';
       gdpeak_signal(i).value(end+1:end+sam) = nan;
    end
end

[val, ~] = max(gdtrough_longd);

for i = 1:size(idx_gd_trough, 2)
    if length(gdtrough_DFFsignal(i).value) < val
       sam = val - length(gdtrough_DFFsignal(i).value); 
       gdtrough_timeDFF(i).time(end+1:end+sam) = length(gdtrough_timeDFF(i).time) + [1:1:sam].';
       gdtrough_DFFsignal(i).value(end+1:end+sam) = nan;
       gdtrough_timestamps(i).time(end+1:end+sam) = length(gdtrough_timestamps(i).time) + [1:1:sam].';
       gdtrough_signal(i).value(end+1:end+sam) = nan;
    end
end

% Transformar structs en matrices (ahora sí con misma longitud):

% PYRAMIDAL:
% Peaks
pyrpeaks = horzcat(pyrpeak_signal.value).';
dfpyrpeaks = horzcat(pyrpeak_DFFsignal.value).';
pyrtimepeaks = [0:size(pyrpeaks, 2)-1]./pyrfs;

figure(21)
subplot(2, 1, 1)
imagesc(pyrtimepeaks, [1:length(pyrpeaks)], pyrpeaks)
subplot(2, 1, 2)
imagesc(pyrtimepeaks, [1:length(pyrpeaks)], dfpyrpeaks)

% Troughs:
pyrtroughs = horzcat(pyrtrough_signal.value).';
dfpyrtroughs = horzcat(pyrtrough_DFFsignal.value).';
pyrtimetroughs = [0:size(pyrtroughs, 2)-1]./pyrfs;

figure(22)
subplot(2, 1, 1)
imagesc(pyrtimetroughs, [1:length(pyrtroughs)], pyrtroughs)
subplot(2, 1, 2)
imagesc(pyrtimetroughs, [1:length(pyrtroughs)], dfpyrtroughs)

% GD:
% Peaks:
gdpeaks = horzcat(gdpeak_signal.value).';
dfgdpeaks = horzcat(gdpeak_DFFsignal.value).';
gdtimepeaks = [0:size(gdpeaks, 2)-1]./gdfs;

figure(23)
subplot(2, 1, 1)
imagesc(gdtimepeaks, [1:length(gdpeaks)], gdpeaks)
subplot(2, 1, 2)
imagesc(gdtimepeaks, [1:length(gdpeaks)], dfgdpeaks)


% Troughs:
gdtroughs = horzcat(gdtrough_signal.value).';
dfgdtroughs = horzcat(gdtrough_DFFsignal.value).';
gdtimetroughs = [0:size(gdtroughs, 2)-1]./gdfs;


figure(24)
subplot(2, 1, 1)
imagesc(gdtimetroughs, [1:length(gdtroughs)], gdtroughs)
subplot(2, 1, 2)
imagesc(gdtimetroughs, [1:length(gdtroughs)], dfgdtroughs)


% Extraer medias y alinear:

% PYR
% Peaks

pyrpeaksmean = nanmean(pyrpeaks, 1); pyrpeaksstd = nanstd(pyrpeaks, 1); 
pyrtroughsmean = nanmean(pyrtroughs,1); pyrtroughsstd = nanstd(pyrtroughs,1);
dfpyrpeaksmean = nanmean(dfpyrpeaks,1); dfpyrpeaksstd = nanstd(dfpyrpeaks,1);
dfpyrtroughsmean = nanmean(dfpyrtroughs,1); dfpyrtroughsstd = nanstd(dfpyrtroughs,1);

gdpeaksmean = nanmean(gdpeaks,1); gdpeaksstd = nanstd(gdpeaks,1); 
gdtroughsmean = nanmean(gdtroughs,1); gdtroughsstd = nanstd(gdtroughs,1);
dfgdpeaksmean = nanmean(dfgdpeaks,1); dfgdpeaksstd = nanstd(dfgdpeaks,1);
dfgdtroughsmean = nanmean(dfgdtroughs,1); dfgdtroughsstd = nanstd(dfgdtroughs,1);

% Pyramidal Peaks:
tiledlayout(2, 1)
nexttile
plot(pyrtimepeaks, pyrpeaksmean ,'m', 'LineWidth', 2)
hold on
plot(pyrtimepeaks, + pyrpeaksstd, 'm', 'LineWidth', 0.3, 'LineStyle', '-')
plot(pyrtimepeaks, - pyrpeaksstd, 'm', 'LineWidth', 0.3,  'LineStyle', '-')
hold off
nexttile
plot(pyrtimepeaks, dfpyrpeaksmean, 'blue', 'LineWidth', 2)
hold on
plot(pyrtimepeaks, + dfpyrpeaksstd, 'b', 'LineWidth', 0.3,  'LineStyle', '-')
plot(pyrtimepeaks, - dfpyrpeaksstd, 'b', 'LineWidth', 0.3,  'LineStyle', '-')
hold off

% Pyramidal troughs:
tiledlayout(2, 1)
nexttile
plot(pyrtimetroughs, pyrtroughsmean ,'m', 'LineWidth', 2)
hold on
plot(pyrtimetroughs, + pyrtroughsstd, 'm', 'LineWidth', 0.3, 'LineStyle', '-')
plot(pyrtimetroughs, - pyrtroughsstd, 'm', 'LineWidth', 0.3,  'LineStyle', '-')
hold off
nexttile
plot(pyrtimetroughs, dfpyrtroughsmean, 'blue', 'LineWidth', 2)
hold on
plot(pyrtimetroughs, + dfpyrtroughsstd, 'b', 'LineWidth', 0.3,  'LineStyle', '-')
plot(pyrtimetroughs, - dfpyrtroughsstd, 'b', 'LineWidth', 0.3,  'LineStyle', '-')
hold off

% GD peaks

tiledlayout(2, 1)
nexttile
plot(gdtimepeaks, gdpeaksmean ,'m', 'LineWidth', 2)
hold on
plot(gdtimepeaks, + gdpeaksstd, 'm', 'LineWidth', 0.3, 'LineStyle', '-')
plot(gdtimepeaks, - gdpeaksstd, 'm', 'LineWidth', 0.3,  'LineStyle', '-')
hold off
nexttile
plot(gdtimepeaks, dfgdpeaksmean, 'blue', 'LineWidth', 2)
hold on
plot(gdtimepeaks, + dfgdpeaksstd, 'b', 'LineWidth', 0.3,  'LineStyle', '-')
plot(gdtimepeaks, - dfgdpeaksstd, 'b', 'LineWidth', 0.3,  'LineStyle', '-')
hold off

% GD troughs

tiledlayout(2, 1)
nexttile
plot(gdtimetroughs, gdtroughsmean ,'m', 'LineWidth', 2)
hold on
plot(gdtimetroughs, + gdtroughsstd, 'm', 'LineWidth', 0.3, 'LineStyle', '-')
plot(gdtimetroughs, - gdtroughsstd, 'm', 'LineWidth', 0.3,  'LineStyle', '-')
hold off
nexttile
plot(gdtimetroughs, dfgdtroughsmean, 'blue', 'LineWidth', 2)
hold on
plot(gdtimetroughs, + dfgdtroughsstd, 'b', 'LineWidth', 0.3,  'LineStyle', '-')
plot(gdtimetroughs, - dfgdtroughsstd, 'b', 'LineWidth', 0.3,  'LineStyle', '-')
hold off


% Cross-correlation:

% Esto es con la amplitud..

tiledlayout(2, 2)
nexttile
crosscorr(double(pyrpeaksmean), double(dfpyrpeaksmean))
ylabel('Pyramidal peaks')
nexttile
crosscorr(double(pyrtroughsmean), double(dfpyrtroughsmean))
ylabel('Pyramidal troughs')
nexttile
crosscorr(double(gdpeaksmean), double(dfgdpeaksmean))
ylabel('GD peaks')
nexttile
crosscorr(double(gdtroughsmean), double(dfgdtroughsmean))
ylabel('GD troughs')



tiledlayout(2, 2)
nexttile
crosscorr(double(pyrpeaksmean(1:100)), double(dfpyrpeaksmean(1:100)))
ylabel('Pyramidal peaks')

nexttile
crosscorr(double(pyrtroughsmean(1:100)), double(dfpyrtroughsmean(1:100)))
ylabel('Pyramidal troughs')

nexttile
crosscorr(double(gdpeaksmean(1:100)), double(dfgdpeaksmean(1:100)))
ylabel('GD peaks')

nexttile
crosscorr(double(gdtroughsmean(1:100)), double(dfgdtroughsmean(1:100)))
ylabel('GD troughs')


% Para hacerlo con el power... downsamplear DFF? 

n = floor(length(DFFZ)/size(pspyr, 2));
dsDF = downsample(DFFZ, n);
if length(dsDF) > size(pspyr, 2)
    dsDF = dsDF(1:size(pspyr, 2));
end

dsDF = double(dsDF);

tiledlayout(3, 1)
nexttile
crosscorr(double(pyrch.meantheta), dsDF)
nexttile
crosscorr(double(pyrch.meangammalow), dsDF)
nexttile
crosscorr(double(pyrch.meangammahigh), dsDF)


n = floor(length(DFFZ)/size(psgd, 2));
dsDFgd = downsample(DFFZ, n);
if length(dsDFgd) > size(psgd, 2)
    dsDFgd = dsDF(1:size(psgd, 2));
end


tiledlayout(3, 1)
nexttile
crosscorr(double(gdch.meantheta), dsDFgd)
nexttile
crosscorr(double(gdch.meangammalow), dsDFgd)
nexttile
crosscorr(double(gdch.meangammahigh), dsDFgd)






%% Guardar todas las figuras:

figpath = strcat(results.FP.path, '/Figures');
figpath1 = strcat(figpath, '/pyrwthresh.jpeg');
figpath2 = strcat(figpath, '/gdwthresh.jpeg');
figpath3 = strcat(figpath, '/pyrposwthresh.jpeg');
figpath4 = strcat(figpath, '/gdposwthresh.jpeg');
figpath5 = strcat(figpath, '/pyrposnothresh.jpeg');
figpath6 = strcat(figpath, '/gdposnothresh.jpeg');
figpath7 = strcat(figpath, '/pyrLFP_FF_wheel.jpeg');
figpath8 = strcat(figpath, '/gdLFP_FF_wheel.jpeg');
figpath9 = strcat(figpath, '/pyrLFP_FF_wheel_smoothed.jpeg');
figpath10 = strcat(figpath, '/gdLFP_FF_wheel_smoothed.jpeg');
figpath11= strcat(figpath, '/pyr_spect.jpeg');
figpath12 = strcat(figpath, '/gd_spect.jpeg');
figpath13 = strcat(figpath, '/pyr_posfq.jpeg');
figpath14 = strcat(figpath, '/gd_posfq.jpeg');
figpath15 = strcat(figpath, '/pyr_posfq_spect.jpeg');
figpath16 = strcat(figpath, '/gd_posfq_spect.jpeg');
figpath17 = strcat(figpath, '/pyrAUC.jpeg');
figpath18 = strcat(figpath, '/gdAUC.jpeg');
figpath19 = strcat(figpath, '/pyrdens.jpeg');
figpath20 = strcat(figpath, '/gddens.jpeg');

saveas(fig1, figpath1);
saveas(fig2, figpath2);
saveas(fig3, figpath3);
saveas(fig4, figpath4);
saveas(fig5, figpath5);
saveas(fig6, figpath6);
saveas(fig7, figpath7);
saveas(fig8, figpath8);
saveas(fig9, figpath9);
saveas(fig10, figpath10);
saveas(fig11, figpath11);
saveas(fig12, figpath12);
saveas(fig13, figpath13);
saveas(fig14, figpath14);
saveas(fig15, figpath15);
saveas(fig16, figpath16);
saveas(fig17, figpath17);
saveas(fig18, figpath18);
saveas(fig19, figpath19);
saveas(fig20, figpath20);




end




