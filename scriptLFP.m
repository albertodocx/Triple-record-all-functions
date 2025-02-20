
% Pipeline for LFP:

% 1. Load and save all files

loadLFPfiles()

% 1.1. Apply 50 Hz notch filter (to process gamma later): 

notch50filter()

% 2. Select channels of interest:

idLFPch()


% Prueba con función de Alberto, para ver si hace lo mismo que lo que he
% hecho: 
% 
% data = results_rawLFP.LFP.LFPdata;
% fs = results_rawLFP.LFP.fs;
% fsdown = fs./1000;
% for i = 1:size(data, 1)
%     newdata(i, :) = downsample(data(i, :), fsdown);
% end
% 
% LCN_compute_psProfile(newdata, 1000, 16)

% 3. Detect theta event:

theta_event()






%% Pruebas load data: 



a = readtable('D:\Users\PC\Desktop\smprueba\Datacontinuous\continuous.dat', 'BinaryType','int16');

mypath = 'C:\Users\PC\OneDrive - Universidad Complutense de Madrid (UCM)\Doctorado\Experimentos\Electro In Vivo\Datos Menopausia\2023-02-28_09-12-26\Record Node 103\experiment1\recording1';
filepath = strcat(mypath, '\structure.oebin');
L = list_open_ephys_binary(filepath,"continuous");
da = load_open_ephys_binary(filepath, "continuous", 1);
B = da.Data(1, 1:10000);

pruebadatos = da.Data(1:16, :).';


    ch = size(LFPdata, 2);
    
    lendata = zeros(1, size(LFPdata, 1));
%     fsdown = fs./1000;
    fsdown = fs./2000; % factor
%     newfs = 1000;
    newfs = 2000;
    lendata = length(downsample(lendata, fsdown));
    dsdata = zeros(lendata, ch);

%     factor_ds = 10;
    disp('downsampling...'); dsdata = downsample(LFPdata,fsdown);

%%

% 3. Extract power spectrum (line plot) of selected channels:


% Esta info ya está contenida en el psprofile que te saca previamente, sólo
% hay que sacar las vars y guardarlas.

ripplech_meanspec = mean(psProfile.(ripple_chname).pS, 2);
thetach_meanspec = mean(psProfile.(theta_chname).pS, 2); % No confundir canal de theta con potencia de theta.
% aquí estoy sacando la media de todas las fq en el canal en el que la
% potencia de theta es mayor (ch16)
fq_spec = psProfile.(ripple_chname).pSfreqs;

idx = fq_spec >= 0 & fq_spec <= 250;

plot(fq_spec(idx), ripplech_meanspec(idx))
hold on
plot(fq_spec(idx), thetach_meanspec(idx))
hold off
legend('ripple channel', 'theta channel')

% 3.1. Interpolate 50 Hz data (power supply noise)




% 4. Theta-delta ratio to detect movement:

% 4.1. Fq thresholds based on Meryl Malezieux, Ashley L. Kees, Christophe Mulle
% (2020)

deltafq = [0.5 3.5];
thetafq = [4 10];

deltaidx = fq_spec >= deltafq(1) & fq_spec <= deltafq(2);
thetaidx = fq_spec >= thetafq(1) & fq_spec <= thetafq(2);

powerfs = 2;
time_spec = [0:size(psProfile.ch5.pS, 2)-1]./powerfs; % 2 porq es la fq de muestreo del pspectrum (2 puntos por cada muestra)

% We look for theta events in the channel where theta power is higher

thetach_ps = psProfile.(theta_chname).pS;

delta = thetach_ps(deltaidx, :);
theta = thetach_ps(thetaidx, :);

delta_mean = mean(delta, 1);
theta_mean = mean(theta, 1);

dt_ratio = delta_mean./theta_mean;

% 4.2. Smooth to find events (absolute values not used later), only indexes are important

dt_ratiosm = movmean(dt_ratio, 20); 

tiledlayout(2, 1)
nexttile
plot(time_spec, dt_ratio)
yline(mean(dt_ratio))
yline(2.5*std(dt_ratio))
nexttile
plot(time_spec, dt_ratiosm)
yline(mean(dt_ratiosm))
yline(2.5*std(dt_ratiosm))


% dt_ratioZ = zscore(dt_ratiosm);
% plot(time_spec, dt_ratioZ)


figure(2)
imagesc(time_spec, fq_spec(thetaidx), theta)
title_spec = sprintf('Mean Power of Theta (%i - %i Hz)', thetafq(1), thetafq(2));
title(title_spec)
colorbar
set(gca, 'YDir', 'normal')


[~,locs,~,~] = findpeaks(dt_ratiosm, powerfs, 'MinPeakProminence', 2*std(dt_ratiosm), 'MinPeakWidth', 2);

plot(time_spec, dt_ratiosm)
yline(2.5*std(dt_ratiosm))
hold on
scatter(locs, repelem(10, length(locs)), '*')

movidx = dt_ratiosm > 2*std(dt_ratiosm);
movevent_theta = nan([1 length(movidx)]);
movevent_theta(movidx) = theta_mean(movidx);
% Mantain all subfq of theta
movevent_pstheta = nan(size(theta));
movevent_pstheta(:, movidx) = theta(:, movidx);

figure(3)
tiledlayout(2, 1)
nexttile
imagesc(time_spec, fq_spec(thetaidx), movevent_pstheta)
title_spec = sprintf('Mean Power of Theta (%i - %i Hz)', thetafq(1), thetafq(2));
title(title_spec)
colorbar
set(gca, 'YDir', 'normal')
nexttile
plot(time_spec, movevent_theta)

% 4.3. Duration of movement event

% Criteria:
% > At least 2 s duration
% > If two events are too closed together (< 1 s) then they are interpreted
% as part of the same event:

movidx = find(movidx ~= 0); % Encuentra valores distintos de 0, donde haya movimiento. 
nonadjacent = [0,find(diff(movidx)>1),length(movidx)]; 
% Al sacar la diff de estos índices me va a decir si la distancia es de 1
% muestra o más. Si es de más muestras, entonces no son
% adyacentes/contiguos. Marca el límite. 
% find(diff(a) > 1) >>> si son mayor que 1, entonces no son contiguos. Genera un vector con los índices en los que se corta (donde no es contiguo). 
% el último valor es length(a)

for k=2:length(nonadjacent)
    bouts{k-1} = movidx(nonadjacent(k-1)+1):movidx(nonadjacent(k)); % Primer valor de los índices que son movimiento hasta el segundo valor en el que deja de ser contiguo (fin del intervalo)
    mov_onset(k-1) = min(bouts{k-1}); % Timestamp onset
    mov_offset(k-1) = max(bouts{k-1}); % Timestamp offset
end

movement_bouts = [mov_onset; mov_offset].';

for i = 1:size(movement_bouts, 1)-1
    if (movement_bouts(i+1, 1) - movement_bouts(i, 2)) <= powerfs % Si es menor o igual a 1 seg (2 muestras para 2 Hz de fs)
       movevents(i, :) = [movement_bouts(i, 1) movement_bouts(i+1, 2)]; 
    else
       movevents(i, :) = [movement_bouts(i, 1) movement_bouts(i, 2)];
    end
end

diffdur = movevents(:, 2) - movevents(:, 1);
threshold = ~(diffdur < powerfs*2); % Menor que 2 seg no interpretamos que sea un evento de movimiento
movevents = movevents(threshold, :);


% 4.4. Extract mean power of theta first during each event and then for the
% whole session: 

for i = 1:size(movevents, 1)
    temp = theta(:, movevents(i, 1):movevents(i, 2));
    temp2 = mean(temp, 1);
    thetaperevent(i) = mean(temp2);
end

thetapersession = mean(thetaperevent);

% Comprobación de valores:

% Con el siguiente código ya me salen los valores bien (viene de código de
% MatLab directamente):
% Eje x: ahora sí es fq en Hz. Eje y: ahora sí que tiene sentido el valor
% de la potencia (no está en el rango de los miles, sino en escala de 0 a
% 50)

% Canal de ripples downsampleado:

rch = results.channels.ch6.dataLFP.';
fs = 1000;
rch_time = ([0:length(rch)-1]./fs).';
rch_time = seconds(rch_time);
data = timetable(rch_time, rch);

[pxx, f] = pspectrum(data);

figure(1)
plot(f,pow2db(pxx))
grid on
xlabel('Frequency (Hz)')
ylabel('Power Spectrum (dB)')
title('Ripple Channel - Default Frequency Resolution')
ylim([0 60])


% Canal de theta downsmapleado:

tch = results.channels.ch16.dataLFP.';
fs = 1000;
tch_time = ([0:length(tch)-1]./fs).';
tch_time = seconds(tch_time);
data = timetable(tch_time, tch);

[pxx2, f] = pspectrum(data);

figure(2)
plot(f,pow2db(pxx2))
grid on
xlabel('Frequency (Hz)')
ylabel('Power Spectrum (dB)')
title('Theta Channel - Default Frequency Resolution')
ylim([0 60])

prueba = pow2db(thetach_ps);
imagesc(time_spec, fq_spec(thetaidx), thetach_ps(thetaidx, :))


%% Comprobaciones: 
% 
% A = psProfile.ch7.pS;
% A = mean(A, 2);
% 
% plot(pxx, A)
% 
% figure(1)
% 
% plot(f, pow2db(pxx), 'blue')
% hold on 
% plot(f, pow2db(A), 'red')
% hold off
% 
% B = psProfile.ch16.pS;
% B = mean(B, 2);
% 
% 
% plot(f, pow2db(pxx2), 'blue')
% hold on 
% plot(f, pow2db(B), 'red')
% hold off

% Sistemáticamente veo que la estimación de la potencia con pspectrum está
% por encima que la estimada con la fx del grupo de Liset (LCN_LFP....etc)

% De momento, me quedo con esa, que es la que está utilizando tb Alberto,
% para que sea igual a como lo calculan ellos. 


%% SACAR RESUMEN DE POWER DE CADA ANIMAL: 

uiwait(msgbox('Select main path (main folder)', 'Instructions', "modal")); 

folderpath = uigetdir(); 


% Get list of all subfolders.
allSubFolders = genpath(folderpath);

% Parse into a cell array.
remain = allSubFolders;
listOfFolderNames = {};
while true
	[singleSubFolder, remain] = strtok(remain, ';');
	if isempty(singleSubFolder)
		break;
	end
	listOfFolderNames = [listOfFolderNames singleSubFolder];
end

lfp_paths = [];
for i = 1:length(listOfFolderNames)
%     if contains(listOfFolderNames(i), 'Data') && contains(listOfFolderNames(i), 'continuous')
    if contains(listOfFolderNames(i), 'recording') && ~contains(listOfFolderNames(i), 'events') && ~contains(listOfFolderNames(i), 'continuous')
        lfp_paths = [lfp_paths listOfFolderNames(i)];
    else
        continue
    end
end


id = '';
expcond = '';
prepost = '';
path = '';

measurenames = {'id', 'expcond', 'prepost', 'path'};

Data = {id, expcond, prepost, path}
