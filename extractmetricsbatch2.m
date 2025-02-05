%%% FUNCIÓN PARA GENERAR MATRIZ DE DATOS DE: 

%%%% RELACIONADAS CON EL TIEMPO %%%%%

% Duración total del freezing 
% Velocidad media del freezing
% Para tiempo medio del episodio de freezing 
% Num de episodios 
% Duración individual de cda episodio de freezing 

% En las siguientes medidas:

% 1. Primer minuto
% 2. Segundo minuto
% 3. Etc
% 4. 
% 5. 

% 6. De los tres primeros minutos 
% 7. De los cinco primeros minutos 

%%%% RELACIONADAS CON LA LOCALIZACIÓN %%%%%

% Velocidad media por zona
% Duración total en cada zona
% Distancia total recorrida por zona
% Episodios totales por zonas
% Duración freezing media en cada zona 
% Distancia mínima al borde (variable continua, en cada t)



% OBJETIVO: Generar *.mat y que permita ir incorporando datos de todos los animales
% (cada fila un animal). El orden de generación de los datos lo elige el
% usuario. 


%%% Extractmetrics.mat (Versión 19.12.2022)

%clear all; clc


function [] = extractmetricsbatch2()

%%%% 

uiwait(msgbox('Select path (a folder/directorio) where you want to save your complete data', 'Instructions', "modal")); 

metricspath = uigetdir(); % Se utiliza esta ruta para guardar todos los *.mat. 


% Prueba para coger muchos archivos: 

uiwait(msgbox('Select main path (main folder)', 'Instructions', "modal"));

mainpath = uigetdir(); 


%%%% ESTO ES DE IMAGEANALYST MATHWORKS:


% Get list of all subfolders.
allSubFolders = genpath(mainpath);

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


videopaths = [];
for i = 1:length(listOfFolderNames)
    if contains(listOfFolderNames(i), 'analyzed')
        videopaths = [videopaths listOfFolderNames(i)];
    else
        continue
    end
end



% Generate seven empty files where information will be stored (arrays):

%%% FILES:

firstdur = 0;
seconddur = 0;
thirddur = 0;
fourthdur = 0;
fifthmin = 0;
velmax = [];

threemin = 0;
fivemin = 0;
path = '';

Data = {firstdur, seconddur, thirddur, fourthdur, fifthmin, threemin, fivemin, velmax, path};

%colnames = {'First Min', 'Second Min', 'Third Min', 'Fourth Min', 'Fifth Min', 'Total 3 Min', 'Total 5 Min', 'Path'}; % Para vars temporales del principio

colnames1 = {'First Min (s)', 'Second Min (s)', 'Third Min (s)', 'Fourth Min (s)', 'Fifth Min (s)', 'Total 3 Min (s) ', 'Total 5 Min (s)', 'Darting', 'Path'}; % Para vars temporales del principio
colnames2 = {'First Min (s)', 'Second Min (s)', 'Third Min (s)', 'Fourth Min (s)', 'Fifth Min (s)', 'Total 3 Min (s) ', 'Total 5 Min (s)', 'Darting', 'Path'}; % Para vars temporales del principio
colnames3 = {'First Min (cm/s)', 'Second Min (cm/s)', 'Third Min (cm/s)', 'Fourth Min (cm/s)', 'Fifth Min (cm/s)', 'Total 3 Min (cm/s) ', 'Total 5 Min (cm/s)', 'Darting', 'Path'}; % Para vars temporales del principio
colnames4 = {'First Min (n)', 'Second Min (n)', 'Third Min (n)', 'Fourth Min (n)', 'Fifth Min (n)', 'Total 3 Min (n) ', 'Total 5 Min (n)', 'Darting', 'Path'}; % Para vars temporales del principio


dmbnames = ["Bin Number", "Time in Bin", "Number Of Episodes", "Total Duration", "Mean Speed in Bin"]; % Para vars relacionadas con distancia mínima al borde + plots (DMB/MDB)
% Esta tabla se puede crear con todas las variables posibles, otra cosa es
% para el plot, si se quieren mostrar unas u otras variables en concreto,
% pero los datos los saca de "todas" (las que digan) las vars de interés
% (de momento DMB x num episodios, DMB x dur total episodios en ese bin,
% DMB x vel media, DMB x tiempo total que pasa en ese bin). 

Data1 = Data;
Data2 = Data;
Data3 = Data;
Data4 = Data;
Data5 = {};
Data6 = {};

filename1 = strcat(metricspath, '/TotalDuration.mat');
filename2 = strcat(metricspath, '/MeanDuration.mat');
filename3 = strcat(metricspath, '/MeanVelocity.mat');
filename4 = strcat(metricspath, '/Episodes.mat');
filename41 = strcat(metricspath, '/TravelledDistance.mat');
%filename5 = strcat(metricspath, '/LocationRelatedMeasures.mat');
filename51 = strcat(metricspath, '/LRM_Speed.mat');
filename52 = strcat(metricspath, '/LRM_Distance.mat');
filename53 = strcat(metricspath, '/LRM_Tottime.mat');
filename54 = strcat(metricspath, '/LRM_Totfrdur.mat');
filename55 = strcat(metricspath, '/LRM_Mboutsdur.mat');
filename56 = strcat(metricspath, '/LRM_Totepi.mat');
filename6 = strcat(metricspath, '/EachEpisodeDuration.mat');
filename7 = strcat(metricspath, '/MinDistancetoBorder.mat');
filename73 = strcat(metricspath, '/MinDistancetoBorder3min.mat');
filename8 = strcat(metricspath, '/MDBxOtherMeasures.mat');
filename83 = strcat(metricspath, '/MDBxOtherMeasures3min.mat');
filename9 = strcat(metricspath, '/MDBxEachEapisodeDuration.mat');
filename93 = strcat(metricspath, '/MDBxEachEpisodeDuration3min.mat');

for i = 1:length(videopaths)

%%%%%%%%%%%%%%%%%%%%% DATOS RELACIONADOS CON EL TIEMPO %%%%%%%%%%%%%%%%%%%%

% Open *.mat file where all relevant variables have been previously saved
% (from BHDtoCSV.m). 

    subpath = strcat(string(videopaths(i)), '/BehData.mat');

    load(subpath)


% Extract indexes for first minute, second minute, etc:
% (siempre van a ser 5 min) 

    minlength = (60.*(length(BHD.Time))./300); % Asumiendo que el total son 5 min/Se puede generalizar más.. 

    firstmin = (1:minlength); 
    secondmin = round(minlength+1:2*minlength);
    thirdmin = round(2*minlength+1:3*minlength);
    fourthmin = round(3*minlength+1:4*minlength);
    fifthmin = round(4*minlength+1:length(BHD.Time));

% Extract other relevant variables for calc:

    Fs = BHD.Fs;
    freez = BHD.freezlogical;
    vel = BHD.MidbackSpeed;
    originalpath = BHD.path;


    velmax = numel(find(BHD.MidbackSpeed>22.9));

%%%%% 1. DURACIÓN TOTAL DEL FREEZING (EN SEGUNDOS): %%%%%

% 1.1. Primeras cinco columnas:

    firstdur = sum(freez(firstmin))./Fs;
    seconddur = sum(freez(secondmin))./Fs;
    thirddur = sum(freez(thirdmin))./Fs;
    fourthdur = sum(freez(fourthmin))./Fs;
    fifthdur = sum(freez(fifthmin))./Fs;

% 1.2. Tres y cinco minutos:

    threemin = (firstdur + seconddur + thirddur)./3; 

    fivemin = (firstdur + seconddur + thirddur + fourthdur + fifthdur)./5;

% 1.3. Añadir a matriz de datos existente (creada previamente) los datos
% obtenidos:

    Data1(i, :) = {firstdur, seconddur, thirddur, fourthdur, fifthdur, threemin, fivemin, velmax, originalpath};



%%%%% 2. TIEMPO MEDIO DE FREEZING (EN SEGUNDOS): %%%%

% ¿Cuántos eventos de freezing hay por cada min? ¿Cuál es la duración media
% de episodio por min?

Bouts  = BHD.Event.Time;

p = 5;
minutes = [0 900 1800 2700 3600 4500];

episodespermin = {};
durboutpermin = {};
epidur = [];

for t = 1:p
    
    rowidx = (Bouts(:, 1) >= minutes(t) & Bouts(:, 1) < minutes(t+1));
    episodes = sum(rowidx);

    if episodes == 0

    dur = 0;

    elseif episodes == 1
        
    allepidur = Bouts(rowidx, 2)./Fs - Bouts(rowidx, 1)./Fs;
    dur = allepidur;
    %dur = (Bouts(rowidx, 2) - Bouts(rowidx, 1))./Fs;
    epidur = [epidur; allepidur];
   
% Si epi es = 1 que coja la duración total de ese episodio

    elseif episodes > 1
    
    allepidur = Bouts(rowidx, 2)./Fs - Bouts(rowidx, 1)./Fs;
    dur = sum(allepidur)./episodes;
    epidur = [epidur; allepidur];

    end

% Guardar datos:
    
    episodespermin{1, t} = episodes;
    durboutpermin{1, t} = dur; % duración media de episodios de freezing

end

    episodespermin = cell2mat(episodespermin);
    durboutpermin = cell2mat(durboutpermin);

    three_meandurbout = mean(epidur(Bouts(:, 1)<2700));
    five_meandurbout = mean(epidur);

%    three_duration = (durboutpermin(:, 1) +  durboutpermin(:, 2) + durboutpermin(:, 3))./3;
%    five_duration = (durboutpermin(:, 1) +  durboutpermin(:, 2) + durboutpermin(:, 3) +  durboutpermin(:, 4) + durboutpermin(:, 5))./5;

    three_episode = (episodespermin(:, 1) + episodespermin(:, 2) + episodespermin(:, 3)); % no está normalizado, es el núm de episodios totales (Y ESTÁ BIEN HECHO)
    five_episode = (episodespermin(:, 1) + episodespermin(:, 2) + episodespermin(:, 3) + episodespermin(:, 4) + episodespermin(:, 5));

    Data2(i, :) = {durboutpermin(:, 1), durboutpermin(:, 2), durboutpermin(:, 3), durboutpermin(:, 4), durboutpermin(:, 5), three_meandurbout, five_meandurbout, velmax, originalpath};
    Data4(i, :) = {episodespermin(:, 1), episodespermin(:, 2), episodespermin(:, 3), episodespermin(:, 4), episodespermin(:, 5), three_episode, five_episode, velmax, originalpath};



%%%% 3. VELOCIDAD MEDIA DEL FREEZING (EN CM/SEG): %%%%

    velfirst = mean(vel(firstmin));
    velsecond = mean(vel(secondmin));
    velthird = mean(vel(thirdmin));
    velfourth = mean(vel(fourthmin));
    velfifth = mean(vel(fifthmin));


    velthree = (velfirst + velsecond + velthird) ./ 3;
    velfive = (velfirst + velsecond + velthird + velfourth + velfifth) ./ 5;


    Data3(i, :) = {velfirst, velsecond, velthird, velfourth, velfifth, velthree, velfive, velmax, originalpath};

    save(filename1, 'Data1');
    save(filename2, 'Data2');
    save(filename3, 'Data3');
    save(filename4, 'Data4');
    
    
%%%%%%%%%%%%%%%%%%%%%%%% DISTANCIA RECORRIDA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dist_tra = BHD.DistanceTraveled; % distancia acumulada en el tiempo

dist_trafirst = dist_tra(firstmin(end));
dist_trasecond = dist_tra(secondmin(end)) - dist_trafirst; % dist en cada min, lo resto con el min anterior para tener el valor absoluto de distancia recorrida en ese min
dist_trathird = dist_tra(thirdmin(end)) - dist_tra(secondmin(end));
dist_trafourth = dist_tra(fourthmin(end)) - dist_tra(thirdmin(end)) ;
dist_trafifth = dist_tra(fifthmin(end)) - dist_tra(fourthmin(end));

dist_tracum3 = dist_tra(thirdmin(end)); % acumulado de los 3 primeros min
dist_tracum5 = dist_tra(fifthmin(end)); % acumulado de lo 5 primeros min

Data41(i, :) = {dist_trafirst, dist_trasecond, dist_trathird, dist_trafourth, dist_trafifth, dist_tracum3, dist_tracum5, originalpath};
measures41 = {'Trav. Dist First Min (cm)', 'Trav. Dist Second Min (cm)', 'Trav. Dist Third Min (cm)', 'Trav. Dist Fourth Min (cm)', 'Trav. Dist Fifth Min (cm)', 'Cum. Dist Three Min (cm)', 'Cum. Dist Five Min (cm)', 'path'};
save(filename41, 'Data41');

%%%%%%%%%%%%%%%%%%% DATOS RELACIONADOS CON LOCALIZACIÓN %%%%%%%%%%%%%%%%%%%

rightlimx = 525; % Límites definitivos 
leftlimx = 125;
highlimy = 440;
lowlimy = 40;

Location = BHD.Location;
locx = Location(:, 1);
locy = Location(:, 2);
locx(locx > rightlimx) = rightlimx;
locx(locx < leftlimx) = leftlimx;
locy(locy > highlimy) = highlimy;
locy(locy < lowlimy) = lowlimy;
Location = [locx locy];


% Establecer límites de centro y de periferia

xaxis_center = [285 365]; % últimos límites acordados con Pablo (19.12.2022)
yaxis_center= [200 280];

xaxis_peri = [125 525]; 
yaxis_peri = [40 440]; 

xaxis_inter = [173 477]; % 3cm * 16 px
yaxis_inter = [88 392];


limcenter = (Location(:, 1) >= xaxis_center(1) & Location(:, 1) <= xaxis_center(2)) & (Location(:, 2) >= yaxis_center(1) & Location(:, 2) <= yaxis_center(2)); 
inter = (((Location(:, 1) >= xaxis_inter(1) & Location(:, 1)  <= xaxis_inter(2))) & ((Location(:, 2) >= yaxis_inter(1) & Location(:, 2) <= yaxis_inter(2)))) & ~limcenter;
peri = (((Location(:, 1) >= xaxis_peri(1) & Location(:, 1) <= xaxis_peri(2))) & ((Location(:, 2) >= yaxis_peri(1) & Location(:, 2) <= yaxis_peri(2)))) & (~limcenter & ~inter) ;
% corner = (BHD.Location(:, 1) <= xaxis_peri(1) | BHD.Location(:, 1) >= xaxis_peri(2)) & (BHD.Location(:, 2) <= yaxis_peri(1) | BHD.Location(:, 2) >= yaxis_peri(2));
% Nota 09.01.2022: las esquinas no están bien calculadas

% Para calcular sólo los 3 primeros minutos:

threeminidx = round(Fs*60*3);

speed3 = BHD.MidbackSpeed(1:threeminidx);
limcenter3 = limcenter(1:threeminidx);
inter3 = inter(1:threeminidx);
peri3 = peri(1:threeminidx);

% Sólo LOCALIZACIÓN:

% 1. Velocidad media del animal en función de centro y periferia y total

%%%%% METRICS.MOVEMENT.DATA = MidbackSpeed. Diferencia es que MidbackSpeed
%%%%% ya está corregido para la muestra que le falta. 

center_speed = mean(BHD.MidbackSpeed(limcenter)); % está en cm/seg
peri_speed = mean(BHD.MidbackSpeed(peri));
inter_speed = mean(BHD.MidbackSpeed(inter));
total_speed = mean(BHD.MidbackSpeed); %por si quieren sacar algún índice más
%corner_speed = mean(BHD.MidbackSpeed(corner));


% 1.1. Velocidad media por zonas en los primeros 3 min:

center_speed3 = mean(BHD.MidbackSpeed(limcenter3)); 
peri_speed3 = mean(BHD.MidbackSpeed(peri3));
inter_speed3 = mean(BHD.MidbackSpeed(inter3));
total_speed3 = mean(BHD.MidbackSpeed);
 
% 2. Duración total del animal en centro, periferia y zona intermedia

center_dur = sum(limcenter)./Fs; % En segundos 
peri_dur = sum(peri)./Fs;
inter_dur = sum(inter)./Fs;
total_dur = sum([center_dur,peri_dur,inter_dur], 'omitnan');
%corner_dur = sum(corner)./Fs;

% 2.1. Duración total por zonas de los 3 primeros min:

center_dur3 = sum(limcenter3)./Fs; % En segundos 
peri_dur3 = sum(peri3)./Fs;
inter_dur3 = sum(inter3)./Fs;
total_dur3 = sum([center_dur3,peri_dur3,inter_dur3], 'omitnan');

% 3. Distancia total recorrida en centro, periferia y total 

% distancia = velocidad media por tiempo 

center_dist = center_speed .* center_dur; % está en cm 
peri_dist = peri_speed .* peri_dur;
inter_dist = inter_speed .* inter_dur;
total_dist = sum([center_dist,peri_dist,inter_dist], 'omitnan');
%corner_dist = corner_speed .* corner_dur;

% 3.1. Tres primeros minutos:

center_dist3 = center_speed3 .* center_dur3; % está en cm 
peri_dist3 = peri_speed3 .* peri_dur3;
inter_dist3 = inter_speed3 .* inter_dur3;
total_dist3 = sum([center_dist3,peri_dist3,inter_dist3], 'omitnan');

% Con FREEZING:

% 4. Episodios totales de freezing en centro y periferia y total

% Centro = 1; Periferia = 0; 
location_center = {}; % Vectores con cada posición durante freezing guardados en cada celda
location_peri = {};
location_inter = {};
%location_corner = {};
location_log = []; % Vector lógico que indica cada freezing si ha sido en el centro o en la periferia
%freez_samples = [];
%freez_cornervector = [];
idx = [];

%  Centro = 1, Peri = 2, Inter = 3. 

for k = 1:size(Bouts, 1)
    location_center{k} = limcenter(Bouts(k, 1):Bouts(k, 2));
    location_peri{k} = peri(Bouts(k, 1):Bouts(k, 2));
    location_inter{k} = inter(Bouts(k, 1):Bouts(k, 2));
    %location_corner{k} = corner(Bouts(k, 1):Bouts(k, 2));
    vector = [sum(location_center{k}), sum(location_peri{k}), sum(location_inter{k})];
    [M(k), idx(k)] = max(vector);
    location_log(k) = idx(k);
end

% Para los 3 minutos:


Bouts3 = Bouts(Bouts(:, 1)<threeminidx, :); 


% Centro = 1; Periferia = 0; 
location_center3 = {}; % Vectores con cada posición durante freezing guardados en cada celda
location_peri3 = {};
location_inter3 = {};
%location_corner3 = {};
location_log3 = []; % Vector lógico que indica cada freezing si ha sido en el centro o en la periferia
%freez_samples3 = [];
%freez_cornervector3 = [];
idx3 = [];

%  Centro = 1, Peri = 2, Inter = 3. 
k = 1;
for k = 1:size(Bouts3, 1)
    if k == size(Bouts3, 1)
        if Bouts3(end, 2) > threeminidx
          location_center3{k} = limcenter3(Bouts3(k, 1):threeminidx);
          location_peri3{k} = peri3(Bouts3(k, 1):threeminidx);
          location_inter3{k} = inter3(Bouts3(k, 1):threeminidx);
          vector = [sum(location_center3{k}), sum(location_peri3{k}), sum(location_inter3{k})];
          [N(k), idx3(k)] = max(vector);
          location_log3(k) = idx3(k);
        end
    else
    location_center3{k} = limcenter3(Bouts3(k, 1):Bouts3(k, 2));
    location_peri3{k} = peri3(Bouts3(k, 1):Bouts3(k, 2));
    location_inter3{k} = inter3(Bouts3(k, 1):Bouts3(k, 2));
    %location_corner3{k} = corner3(Bouts3(k, 1):Bouts3(k, 2));
    vector = [sum(location_center3{k}), sum(location_peri3{k}), sum(location_inter3{k})];
    [N(k), idx3(k)] = max(vector);
    location_log3(k) = idx3(k);
    end
end


%% LRM RESPECTO AL FREEZING


% for k = 1:size(Bouts, 1)
%    if sum(location_corner{k}) > length(location_corner{k})./2
%     freez_cornervector(k) = 1;
%    elseif sum(location_corner{k}) <= length(location_corner{k})./2
%     freez_cornervector(k) = 0;
%    end
% end

% Episodios totales por zonas:

totepi_center = sum(idx == 1);
totepi_peri = sum(idx == 2);
totepi_inter = sum(idx == 3);
%totepi_corner = sum(freez_cornervector);
totepi_total = (totepi_center + totepi_peri + totepi_inter);

totepi_centerper = totepi_center/totepi_total*100;
totepi_periper = totepi_peri/totepi_total* 100;
totepi_interper = totepi_inter/totepi_total *100;

% Epis por zonas 3 min:

totepi_center3 = sum(idx3 == 1);
totepi_peri3 = sum(idx3 == 2);
totepi_inter3 = sum(idx3 == 3);
%totepi_corner = sum(freez_cornervector);
totepi_total3 = (totepi_center3 + totepi_peri3 + totepi_inter3);


totepi_centerper3 = totepi_center3/totepi_total3 *100;
totepi_periper3 = totepi_peri3/totepi_total3 * 100;
totepi_interper3 = totepi_inter3/totepi_total3 *100;

% 5. Duración total 5 min  freezing:

boutsdur = (Bouts(:, 2) - Bouts(:, 1) + 1)./Fs; % Para igualarlo a lo que sale del bucle anterior

boutsdur = boutsdur.';

freezdur_center = sum(boutsdur(idx == 1)); % 
freezdur_peri = sum(boutsdur(idx == 2));
freezdur_inter = sum(boutsdur(idx == 3));
freezdur_total = mean(boutsdur); % Duración media del freezing independientemente de la zona. 
%freezdur_corner = sum(boutsdur(freez_cornervector == 1))./ sum(freez_cornervector);

% Duración total 3 min:

freez3 = freez(1:threeminidx);

freezdur_center3 = sum(freez3(freez3 == 1 & limcenter3 == 1))./Fs;
freezdur_peri3 = sum(freez3(freez3 == 1 & inter3 == 1))./Fs;
freezdur_inter3 = sum(freez3(freez3 == 1 & peri3 == 1))./Fs;
freezdur_total3 = (freezdur_center3 + freezdur_peri3 + freezdur_inter3);


% 6. Duración freezing media 

mfreezdur_center = sum(boutsdur(idx == 1))./sum(idx == 1); % Duración media de freezing en el centro >>> duración total/n episodios en esa zona
mfreezdur_peri = sum(boutsdur(idx == 2))./sum(idx == 2);
mfreezdur_inter = sum(boutsdur(idx == 3))./sum(idx == 3);
mfreezdur_total = mean(boutsdur); % Duración media del freezing independientemente de la zona. 
%freezdur_corner = sum(boutsdur(freez_cornervector == 1))./ sum(freez_cornervector);


% Por 3 min

boutsdur3 = (Bouts3(:, 2) - Bouts3(:, 1) + 1)./Fs; % Para igualarlo a lo que sale del bucle anterior

boutsdur3 = boutsdur3.';

mfreezdur_center3 = sum(boutsdur3(idx3 == 1))./sum(idx3 == 1); % Duración media de freezing en el centro >>> duración total/n episodios en esa zona
mfreezdur_peri3 = sum(boutsdur3(idx3 == 2))./sum(idx3 == 2);
mfreezdur_inter3 = sum(boutsdur3(idx3 == 3))./sum(idx3 == 3);
mfreezdur_total3 = mean(boutsdur3); % Duración media del freezing independientemente de la zona. 


% Calculate location measures ratios:

if isnan(center_speed)
    ratio_v = 0;
elseif isnan(peri_speed)
    ratio_v = 1;
else
    ratio_v = center_speed/(center_speed + peri_speed);
end

if isnan(center_dur)
    ratio_t = 0;
elseif isnan(peri_dur)
    ratio_t = 1;
else
    ratio_t = center_dur/(center_dur + peri_dur);
end

if isnan(center_dist)
    ratio_d = 0;
elseif isnan(peri_dist)
    ratio_d = 1;
else
    ratio_d = center_dist/(center_dist + peri_dist);
end

if isnan(totepi_center)
    ratio_e = 0;
elseif isnan(totepi_peri)
    ratio_e = 1;
else
    ratio_e = totepi_center/(totepi_center + totepi_peri);
end

if isnan(freezdur_center)
    ratio_ft = 0;
elseif isnan(freezdur_peri)
    ratio_ft = 1;
else
    ratio_ft = freezdur_center/(freezdur_center + freezdur_peri);
end

if isnan(mfreezdur_center)
    ratio_mft = 0;
elseif isnan(mfreezdur_peri)
    ratio_mft = 1;
else
    ratio_mft = mfreezdur_center/(mfreezdur_center + mfreezdur_peri);
end

% Ratios for 3 min:

if isnan(center_speed3)
    ratio_v3 = 0;
elseif isnan(peri_speed3)
    ratio_v3 = 1;
else
    ratio_v3 = center_speed3/(center_speed3 + peri_speed3);
end

if isnan(center_dur3)
    ratio_t3 = 0;
elseif isnan(peri_dur3)
    ratio_t3 = 1;
else
    ratio_t3 = center_dur3/(center_dur3 + peri_dur3);
end

if isnan(center_dist3)
    ratio_d3 = 0;
elseif isnan(peri_dist3)
    ratio_d3 = 1;
else
    ratio_d3 = center_dist3/(center_dist3 + peri_dist3);
end

if isnan(totepi_center3)
    ratio_e3 = 0;
elseif isnan(totepi_peri3)
    ratio_e3 = 1;
else
    ratio_e3 = totepi_center3/(totepi_center3 + totepi_peri3);
end

if isnan(freezdur_center3)
    ratio_ft3 = 0;
elseif isnan(freezdur_peri3)
    ratio_ft3 = 1;
else
    ratio_ft3 = freezdur_center3/(freezdur_center3 + freezdur_peri3);
end

if isnan(mfreezdur_center3)
    ratio_mft3 = 0;
elseif isnan(mfreezdur_peri3)
    ratio_mft3 = 1;
else
    ratio_mft3 = mfreezdur_center3/(mfreezdur_center3 + mfreezdur_peri3);
end

% GUARDAR NOMBRES DE VARS:

% %measurenames = {'M_Speed Center (cm/s)', 'M_Speed Inter (cm/s)', 'M_Speed Peri (cm/s)', 'M_Speed Total (cm/s)',... 
%  %           'Tot_Dur Center (s)', 'Tot_Dur Inter (s)', 'Tot_Dur Peri (s)', 'Tot_Dur Total (s)',... 
%   %          'Tot_Dist Center (cm)', 'Tot_Dist Inter (cm)', 'Tot_Dist Peri (cm)', 'Tot_Dist Total (cm)',...
% %
%  %           'Tot_Epi Center (n)', 'Tot_Epi Inter (n)', 'Tot_Epi Peri (n)', 'Tot_Epi Total (n)',...
%             'M_FreezDur Center (s)','M_FreezDur Inter (s)','M_FreezDur Peri (s)','M_FreezDur Total (s)',...
%          'V-Ratio', 'T-Ratio','D-Ratio','E-Ratio','FT-Ratio','path'};
% 

measurenames51 = {'M_Speed Center 3min (cm/s)', 'M_Speed Inter 3min (cm/s)', 'M_Speed Peri 3min (cm/s)', 'M_Speed Total 3min (cm/s)', 'V-Ratio 3 min',...
                 'M_Speed Center 5min (cm/s)', 'M_Speed Inter 5min (cm/s)', 'M_Speed Peri 5min (cm/s)', 'M_Speed Total 5min (cm/s)', 'V-Ratio 5 min', 'path'};
measurenames52 = {'Tot_time Center 3min (s)', 'Tot_time Inter 3min (s)', 'Tot_time Peri 3min (s)', 'Tot_time Total 3min (s)', 'T-Ratio 3min',...
                  'Tot_time Center 5min (s)', 'Tot_time Inter 5min (s)', 'Tot_time Peri 5min (s)', 'Tot_time Total 5min (s)', 'T-Ratio 5min', 'path'};
measurenames53 = {'Tot_Dist Center 3min (cm)', 'Tot_Dist Inter 3min (cm)', 'Tot_Dist Peri 3min (cm)', 'Tot_Dist Total 3min (cm)', 'D-Ratio 3min',...
                  'Tot_Dist Center 5min (cm)', 'Tot_Dist Inter 5min (cm)', 'Tot_Dist Peri 5min (cm)', 'Tot_Dist Total 5min (cm)', 'D-Ratio 5min','path'};

measurenames54 = {'Tot_frdur Center 3min (s)','Tot_frdur Inter 3min (s)','Tot_frdur Peri 3min (s)','Tot_frdur Total 3min (s)','FT-Ratio 3min',...
                 'Tot_frdur Center 5min (s)','Tot_frdur Inter 5min (s)','Tot_frdur Peri 5min (s)','Tot_frdur Total 5min (s)','FT-Ratio 5min', 'path'};
measurenames55 = {'M_FreezDur Center 3min (s)','M_FreezDur Inter 3min (s)','M_FreezDur Peri 3min (s)','M_FreezDur Total 3min (s)','MFT-Ratio 3min',...
                  'M_FreezDur Center 5min (s)','M_FreezDur Inter 5min (s)','M_FreezDur Peri 5min (s)','M_FreezDur Total 5min (s)','MFT-Ratio 5min', 'path'};
measurenames56 = {'Tot_Epi Center 3min (n)', 'Tot_Epi Inter 3min (n)', 'Tot_Epi Peri 3min (n)', 'Tot_Epi Total 3min (n)','E-Ratio 3min', '% Tot_Epi Center 3min (%)', '% Tot_Epi Inter 3min (%)', '% Tot_Epi Peri 3min (%)',...
                 'Tot_Epi Center 5min (n)', 'Tot_Epi Inter 5min (n)', 'Tot_Epi Peri 5min (n)', 'Tot_Epi Total 5min (n)','E-Ratio 5min',  '% Tot_Epi Center 5min (%)', '% Tot_Epi Inter 5min (%)', '% Tot_Epi Peri 5min (%)','path'};

% Sin esquina:
% Data5(i, :) = {center_speed, inter_speed, peri_speed, total_speed,...
%                center_dur, inter_dur, peri_dur, total_dur,...
%                 center_dist, inter_dist, peri_dist, total_dist,...
%                 totepi_center, totepi_inter, totepi_peri, totepi_total,...
%                 freezdur_center, freezdur_inter, freezdur_peri, freezdur_total,...
%                 ratio_v, ratio_t, ratio_d, ratio_e, ratio_ft,...
%                 originalpath};


% LAS SEIS VARIABLES DIVIDIDAS EN DISTINTOS EXCELS: velocidad en cada zona,
% distancia en cada zona, duración total en cada zona + duración total del
% freezing en cada zona, duración media de episodios de freezing en cada
% zona y num de episodios por cada zona (medidas: 5 min + ratio + 3 min +
% ratio)

% Velocidad en loc:

Data51(i, :) = {center_speed3, inter_speed3, peri_speed3, total_speed3, ratio_v3,...
               center_speed, inter_speed, peri_speed, total_speed, ratio_v, originalpath};
% Duración en cada loc:
Data52(i, :) = {center_dur3, inter_dur3, peri_dur3, total_dur3, ratio_t3,...
               center_dur, inter_dur, peri_dur, total_dur, ratio_t, originalpath};
% Distancia en cada loc:
Data53(i, :) = {center_dist3, inter_dist3, peri_dist3, total_dist3,ratio_d3,...
                center_dist, inter_dist, peri_dist, total_dist,ratio_d, originalpath};
% Duración total de episodios en cada loc:
Data54(i, :) = {freezdur_center3, freezdur_inter3, freezdur_peri3, freezdur_total3, ratio_ft3,...
                freezdur_center, freezdur_inter, freezdur_peri, freezdur_total,ratio_ft, originalpath};
% Duración media de episodios en cada loc:
Data55(i, :) = {mfreezdur_center3, mfreezdur_inter3, mfreezdur_peri3, mfreezdur_total3, ratio_mft3,...
                mfreezdur_center, mfreezdur_inter, mfreezdur_peri, mfreezdur_total,ratio_mft, originalpath};
% Episodios totales de freezing en cada loc:
Data56(i, :) = {totepi_center3, totepi_inter3, totepi_peri3, totepi_total3,ratio_e3, totepi_centerper3, totepi_interper3, totepi_periper3,...
               totepi_center, totepi_inter, totepi_peri, totepi_total,ratio_e,totepi_centerper, totepi_interper, totepi_periper, originalpath};



Data6(i, :) = {boutsdur}; % Guardo duración de cada episodio de freezing por vídeo

    save(filename51, 'Data51');
    save(filename52, 'Data52');
    save(filename53, 'Data53');
    save(filename54, 'Data54');
    save(filename55, 'Data55');
    save(filename56, 'Data56');
    save(filename6, 'Data6'); 





%%% GUARDAR VARIABLES EN WORKSPACE Y EN MAT ESQUINA

% Esto igual se puede añadir todo en el mismo Excel, pero en distintas
% hojas de Excel. 

% Con la esquina
% measurenames = {'M_Speed Center', 'M_Speed Inter', 'M_Speed Peri', 'M_Speed Total', 'M_Speed Corner',... 
%             'Tot_Dur Center', 'Tot_Dur Inter', 'Tot_Dur Peri', 'Tot_Dur Total','Tot_Dur Corner',... 
%             'Tot_Dist Center', 'Tot_Dist Inter', 'Tot_Dist Peri', 'Tot_Dist Total','Tot_Dist Corner',...
%             'Tot_Epi Center', 'Tot_Epi Inter', 'Tot_Epi Peri', 'Tot_Epi Total','Tot_Epi Corner',...
%             'M_FreezDur Center','M_FreezDur Inter','M_FreezDur Peri','M_FreezDur Total','M_FreezDur Corner','path'};
% 
%  Sin esquina:

% Con la esquina
% Data5(i, :) = {center_speed, inter_speed, peri_speed, total_speed, corner_speed,...
%                center_dur, inter_dur, peri_dur, total_dur, corner_dur,...
%                 center_dist, inter_dist, peri_dist, total_dist, corner_dist,...
%                 totepi_center, totepi_inter, totepi_peri, totepi_total, totepi_corner,...
%                 freezdur_center, freezdur_inter, freezdur_peri, freezdur_total, freezdur_corner,...
%                 originalpath};




% 5. Distancia mínima al borde (el que sea) en cada punto del registro.
% Output: vector de misma longitud que pos con el valor en cm de la
% distancia mínima al borde. 

rightlimx = 525; % Límites definitivos 
leftlimx = 125;
highlimy = 440;
lowlimy = 40;

Location = BHD.Location;
locx = Location(:, 1);
locy = Location(:, 2);
locx(locx > rightlimx) = rightlimx;
locx(locx < leftlimx) = leftlimx;
locy(locy > highlimy) = highlimy;
locy(locy < lowlimy) = lowlimy;
Location = [locx locy];

% Comprobación de que está bien: 
% min(Location(:, 1))
% max(Location(:, 1))
% max(Location(:, 2))
% min(Location(:, 2))


tmp = [];
minVal = [];
minDist = [];
frompxtocm = 16; % 16 pixeles/cm
rightdist = [];
leftdist = [];
highdist = [];
lowdist = [];
for k = 1:length(Location)
    rightdist(k) = rightlimx - Location(k, 1);
    leftdist(k) = Location(k, 1) - leftlimx;
    highdist(k) = highlimy - Location(k, 2);
    lowdist(k) = Location(k, 2) - lowlimy;
    tmp = [rightdist(k), leftdist(k), highdist(k), lowdist(k)];
    minVal = min(tmp);
    minDist(k) = minVal/frompxtocm;
end


maxd = 12.5; % distancia max que es posible tener --- max está bien calculado, ya he comprobado y el max que sale del bucle es 12.5
mind = 0; % distancia mín que es posible tener 
minDist(minDist < 0) = 0;
minDist(minDist > maxd) = maxd;

minDist3 = minDist(1:threeminidx);


%%% >> 5.1. Extraer datos de DMB x varcont:

xlim2 = [0:(maxd*0.1):maxd];
% 5.1.1. Tiempo total que está en ese bin: 

bindist = {};
for k = 1:(length(xlim2)-1)
    bindist{i, k} = find(minDist >= xlim2(k) & minDist < xlim2(k+1)); % Índices de cuando el animal está en cada bin (más cerca o más lejos del borde)
    bindist_time(i, k) = length(bindist{i, k})./Fs; % Tiempo que pasa en cada bin
end

% Para 3 min:


bindist3 = {};
for k = 1:(length(xlim2)-1)
    bindist3{i, k} = find(minDist3 >= xlim2(k) & minDist3 < xlim2(k+1)); % Índices de cuando el animal está en cada bin (más cerca o más lejos del borde)
    bindist_time3(i, k) = length(bindist3{i, k})./Fs; % Tiempo que pasa en cada bin
end




% 5.1.2. Distancia total recorrida en cada bin: 

distancetraveled = BHD.DistanceTraveled;
distancetraveled = diff(distancetraveled);
distancetraveled = [0; distancetraveled];

bindistance = [];
for k = 1:(length(xlim2)-1)
    bindistancecell{i, k} = distancetraveled(minDist >= xlim2(k) & minDist < xlim2(k+1)); % Índices de cuando el animal está en cada bin (más cerca o más lejos del borde)
    bindistance(i, k) = sum(bindistancecell{i, k}); % Distancia total recorrida en cada bin (sumatorio de las diferencias en cada punto de muestreo) 
end

% Para 3 min: 

distancetraveled3 = distancetraveled(1:threeminidx);
bindistance3 = [];
for k = 1:(length(xlim2)-1)
    bindistancecell3{i, k} = distancetraveled3(minDist3 >= xlim2(k) & minDist3 < xlim2(k+1)); % Índices de cuando el animal está en cada bin (más cerca o más lejos del borde)
    bindistance3(i, k) = sum(bindistancecell3{i, k}); % Distancia total recorrida en cada bin (sumatorio de las diferencias en cada punto de muestreo) 
end

% 5.1.3. ¿Cuántos eventos de freezing en cada bin?


% De momento nos quedamos sólo con el onset (habría que revisar en datos
% "fiables" si durante el episodio de freezing realmente se mueven tanto o
% no)

mindistfreez = minDist(Bouts(:, 1)); % Lo utilizo para poder etiquetar después cada duración del freezing con su correspondiente bin
mindistfreez = mindistfreez.';

binfreezepi = [] ; % Este es un vector en el que cada índice representa cada bin del espacio (de pegado al borde a en el centro)
for k = 1:(length(xlim2) - 1)
    binfreezepi(i, k) = sum(mindistfreez >= xlim2(k) & mindistfreez < xlim2(k+1));
    
end

binfreezepinorm(i, :) = binfreezepi(i,:)./sum(binfreezepi(i, :));

% Para 3  min: 

mindistfreez3 = minDist3(Bouts3(:, 1)); % Lo utilizo para poder etiquetar después cada duración del freezing con su correspondiente bin
mindistfreez3 = mindistfreez3.';

binfreezepi3 = [] ; % Este es un vector en el que cada índice representa cada bin del espacio (de pegado al borde a en el centro)
for k = 1:(length(xlim2) - 1)
    binfreezepi3(i, k) = sum(mindistfreez3 >= xlim2(k) & mindistfreez3 < xlim2(k+1));
    
end

binfreezepinorm3(i, :) = binfreezepi3(i,:)./sum(binfreezepi3(i, :));


% 5.1.4. ¿Cuánto dura cada evento de freezing dentro de cada bin? y ¿
% Cuánto dura en total el freezing del animal que ha hecho en un bin? 

% Etiqueta en cada freezing:
% 5 min:
boutsdur = boutsdur.';
for k = 1:length(mindistfreez)
    for j = 1:(length(xlim2) - 1)
        if (mindistfreez(k) >= xlim2(j) && mindistfreez(k) < xlim2(j+1))
            boutsdur(k, 2) = xlim2(j+1); 
        else
            continue
        end
    end
end 

% Para 3 min:

boutsdur3 = boutsdur3.';
for k = 1:length(mindistfreez3)
    for j = 1:(length(xlim2) - 1)
        if (mindistfreez3(k) >= xlim2(j) & mindistfreez3(k) < xlim2(j+1))
            boutsdur3(k, 2) = xlim2(j+1); % EL + 0,5 es simplemente para que se puedan representar en el medio del bin los valores correspondientes (entre 1 y 2 por ej) y para que no aparezcan alineados con el valor en x.... es por motivos puramente prácticos
        else
            continue
        end
    end
end 

% 5 min:

binmeandurepi = zeros([length(videopaths) (length(xlim2)-1)]);
bintotaldurepi = zeros([length(videopaths) (length(xlim2)-1)]);

if isempty(boutsdur)
    message = strcat('There are no freezing episodes for this file: ', videopaths(i));
    warndlg(message)
    
    for k = 1:(length(xlim2) -1)
        bin_name = sprintf('Bin_%i', k);
        catnames(k) = string(bin_name); 
        bineachdurepi(i).(bin_name) = 0;
    end

    bintotaldurepi(i, 1:length(xlim2)-1) = zeros([1 length(xlim2)-1]);
    binmeandurepi(i, 1:length(xlim2)-1) = zeros([1 length(xlim2)-1]);
    
else
    for k = 1:(length(xlim2) - 1)

     % Los que estén clasificados con el mismo núm se agrupan (si es dur total se suma, si es dur individual se guardan en forma de tabla dinámica porq van a tener distintas longitudes los vectores)
        idxtemp = find(boutsdur(:, 2) == xlim2(k + 1));
           
        bintotaldurepi(i, k) = sum(boutsdur(idxtemp, 1)); % Igual aquí hay que cambiarlo para que se quede cada animal en cada fila (bintotaldurepi(i, k))
      
        bin_name = sprintf('Bin_%i', k);
        catnames(k) = string(bin_name); 
        bineachdurepi(i).(bin_name) = boutsdur(idxtemp); % Genera tabla dinámica en la que cada celda puede ser un vector con la duración de cada episodio que ha tenido el animal por bin. 
        
        binmeandurepi(i, k) = mean(bineachdurepi(i).(bin_name));
        
        % Duración total y núm episodios de freezing totales normalizados:

        
    end
end


bintotaldurepinorm(i, :) = bintotaldurepi(i, :)./sum(bintotaldurepi(i, :));

% Para 3 min: 

binmeandurepi3 = zeros([length(videopaths) (length(xlim2)-1)]);
bintotaldurepi3 = zeros([length(videopaths) (length(xlim2)-1)]);

if isempty(boutsdur3)
    message = strcat('There are no freezing episodes for this file: ', videopaths(i));
    warndlg(message)
    
    for k = 1:(length(xlim2) -1)
        bin_name3 = sprintf('Bin_%i', k);
        catnames3(k) = string(bin_name3); 
        bineachdurepi3(i).(bin_name3) = 0;
    end

    bintotaldurepi3(i, 1:length(xlim2)-1) = zeros([1 length(xlim2)-1]);
    binmeandurepi3(i, 1:length(xlim2)-1) = zeros([1 length(xlim2)-1]);
    
else
    for k = 1:(length(xlim2) - 1)

     % Los que estén clasificados con el mismo núm se agrupan (si es dur total se suma, si es dur individual se guardan en forma de tabla dinámica porq van a tener distintas longitudes los vectores)
        idxtemp3 = find(boutsdur3(:, 2) == xlim2(k + 1));
           
        bintotaldurepi3(i, k) = sum(boutsdur3(idxtemp3, 1)); % Igual aquí hay que cambiarlo para que se quede cada animal en cada fila (bintotaldurepi(i, k))
      
        bin_name3 = sprintf('Bin_%i', k);
        catnames3(k) = string(bin_name3); 
        bineachdurepi3(i).(bin_name3) = boutsdur3(idxtemp3); % Genera tabla dinámica en la que cada celda puede ser un vector con la duración de cada episodio que ha tenido el animal por bin. 
        
        binmeandurepi3(i, k) = mean(bineachdurepi3(i).(bin_name3));
        
        % Duración total y núm episodios de freezing totales normalizados:

        
    end
end


bintotaldurepinorm3(i, :) = bintotaldurepi3(i, :)./sum(bintotaldurepi3(i, :));


% 5.1.5. ¿Cuál es la velocidad media en cada bin?

% 5min:

binmeanspeed = [];

for k = 1:(length(xlim2)-1)
    binmeanspeed(i, k) = mean(vel(find(minDist >= xlim2(k) & minDist < xlim2(k+1)))); % Comprobar que hace bien lo que quiero que haga ...
end

% 3min:

vel3 = vel(1:threeminidx);
binmeanspeed3 = [];

for k = 1:(length(xlim2)-1)
    binmeanspeed3(i, k) = mean(vel3(find(minDist3 >= xlim2(k) & minDist3 < xlim2(k+1)))); % Comprobar que hace bien lo que quiero que haga ...
end

Data7(i, :) = {minDist};
save(filename7, 'Data7'); 

Data8(i, :) = {catnames.', bindist_time(i, :).', binmeanspeed(i, :).', bindistance(i, :).', binfreezepi(i, :).', bintotaldurepi(i, :).', binmeandurepi(i, :).', binfreezepinorm(i, :).', bintotaldurepinorm(i, :).'};
% Data8(i, :) = array2table( [bindist_time(i, :).', binmeanspeed(i, :).', bindistance(i, :).', binfreezepi(i, :).', bintotaldurepi(i, :).', binmeandurepi(i, :).', binfreezepinorm(i, :).', bintotaldurepinorm(i, :).'], ...
%     'VariableNames', ["BinDist Time (s)", "BinMean Speed (cm/s)", "BinDistance (cm)", "BinFreezEpisodes (n)", "BinTotalDurEpi", "BinMeanDurEpi (s)", "BinFreezeEpiNorm (s/min)", "BinTotalDurepiNorm (s/min)"]);
% catnames = table(catnames.','VariableNames', "Bin Number");
% Data8 = [catnames];


save(filename8, 'Data8')
Data9(i, :) = {bineachdurepi};
save(filename9, 'Data9')

Data73(i, :) = {minDist3};
save(filename73, 'Data73'); 

Data83(i, :) = {catnames3.', bindist_time3(i, :).', binmeanspeed3(i, :).', bindistance3(i, :).', binfreezepi3(i, :).', bintotaldurepi3(i, :).', binmeandurepi3(i, :).', binfreezepinorm3(i, :).', bintotaldurepinorm3(i, :).'};
save(filename83, 'Data83')
Data93(i, :) = {bineachdurepi3};
save(filename93, 'Data93')


end % Fin del bucle


%  allminDist = cell2mat(Data7);
 
allminDist = [];

  for nn = 1:size(Data7, 1)
      if length(Data7{nn}) < 4500
          Data7{nn} = [Data7{nn}, nan(1, 4500 - length(Data7{nn}))];
          allminDist(nn, :) = Data7{nn};
      elseif length(Data7{nn}) > 4500
          Data7{nn} = Data7{nn}(1:4500);
          allminDist(nn, :) = Data7{nn};
      end
  end
  
  allminDist3 = [];

  for nn = 1:size(Data73, 1)
      if length(Data73{nn}) < 4500
          Data73{nn} = [Data73{nn}, nan(1, 4500 - length(Data73{nn}))];
          allminDist(nn, :) = Data73{nn};
      elseif length(Data73{nn}) > 4500
          Data73{nn} = Data73{nn}(1:4500);
          allminDist(nn, :) = Data73{nn};
      end
  end

% Guardar las tablas en formato XLSX para que se pueda leer en excel:

    Data1 = array2table(Data1, 'VariableNames', colnames1);
    Data2 = array2table(Data2, 'VariableNames', colnames2);
    Data3 = array2table(Data3, 'VariableNames', colnames3);
    Data4 = array2table(Data4, 'VariableNames', colnames4);
    Data41 = array2table(Data41, 'VariableNames', measures41);
%     Data5 = array2table(Data5, 'VariableNames', measurenames);
    Data51 = array2table(Data51, 'VariableNames', measurenames51);
    Data52 = array2table(Data52, 'VariableNames', measurenames52);
    Data53 = array2table(Data53, 'VariableNames', measurenames53);
    Data54 = array2table(Data54, 'VariableNames', measurenames54);
    Data55 = array2table(Data55, 'VariableNames', measurenames55);
    Data56 = array2table(Data56, 'VariableNames', measurenames56);
    Data6 = array2table(Data6);
     Data7 = array2table(Data7);
     Data73 = array2table(Data73);
%      Data8 = array2table(Data8);
%      Data83 = array2table(Data83);
     Data9 = array2table(Data9);
     Data93 = array2table(Data93);


    xlsx1 = strcat(metricspath, '/TotalDuration.xlsx');
    xlsx2 = strcat(metricspath, '/MeanDuration.xlsx');
    xlsx3 = strcat(metricspath, '/MeanVelocity.xlsx');
    xlsx4 = strcat(metricspath, '/Episodes.xlsx');
    xlsx41 = strcat(metricspath, '/TravelledDistance.xlsx');
    %xlsx5 = strcat(metricspath, '/LocationRelatedMeasures.xlsx');
    xlsx51 = strcat(metricspath, '/LRM_Speed.xlsx');
    xlsx52 = strcat(metricspath, '/LRM_TotTime.xlsx');
    xlsx53 = strcat(metricspath, '/LRM_TotDist.xlsx');
    xlsx54 = strcat(metricspath, '/LRM_Totfrdur.xlsx');
    xlsx55 = strcat(metricspath, '/LRM_Mboutsdur.xlsx');
    xlsx56 = strcat(metricspath, '/LRM_Totepi.xlsx');
    xlsx6 = strcat(metricspath, '/EachEpisodeDuration.xlsx');
     xlsx7 = strcat(metricspath, '/MinDistancetoBorder.xlsx');
     xlsx73 = strcat(metricspath, '/MinDistancetoBorder3min.xlsx');
     xlsx8 = strcat(metricspath, '/MDBandothervars.xlsx');
     xlsx83 = strcat(metricspath, 'MDBandothervars3min.xlsx');
     xlsx9 = strcat(metricspath, '/MDBbineachepidur.xlsx');
     xlsx93 = strcat(metricspath, '/MDBbineachepidur3min.xlsx');

    writetable(Data1, xlsx1, 'WriteVariableNames', true); 
    writetable(Data2, xlsx2, 'WriteVariableNames', true); 
    writetable(Data3, xlsx3, 'WriteVariableNames', true); 
    writetable(Data4, xlsx4, 'WriteVariableNames', true); 
    writetable(Data41, xlsx41, 'WriteVariableNames', true);
    %writetable(Data5, xlsx5, 'WriteVariableNames', true); 
    writetable(Data51, xlsx51, 'WriteVariableNames', true); 
    writetable(Data52, xlsx52, 'WriteVariableNames', true); 
    writetable(Data53, xlsx53, 'WriteVariableNames', true); 
    writetable(Data54, xlsx54, 'WriteVariableNames', true); 
    writetable(Data55, xlsx55, 'WriteVariableNames', true); 
    writetable(Data56, xlsx56, 'WriteVariableNames', true); 
    writetable(Data6, xlsx6, 'WriteVariableNames', 0);
      writematrix(allminDist, xlsx7);
      writematrix(allminDist3, xlsx73);
     writetable(Data9, xlsx9, 'WriteVariableNames', 0);
     writetable(Data93, xlsx93, 'WriteVariableNames', 0);
    

% Generar y guardar plots de distancia mínima al borde y todas las vars que
% se puedan/quieran asociar


%%%%%%%%         REPRESENTACIÓN GRÁFICA             %%%%%%%%%%%%%%%%%%%

% Reconfiguro datos para poder trabajar con boxchart:
% 

bin5minquant = cell2mat(Data8(:, 2:end)); bin5minquant = array2table(bin5minquant, 'VariableNames', ["BinDist Time (s)", "BinMean Speed (cm/s)", "BinDistance (cm)", "BinFreezEpisodes (n)", "BinTotalDurEpi", "BinMeanDurEpi (s)", "BinFreezeEpiNorm (s/min)", "BinTotalDurepiNorm (s/min)"]);
bin5minstr = repmat(catnames.', [length(videopaths) 1]); bin5minstr = array2table(bin5minstr, 'VariableNames', "Bin Number");
mdbdata = [bin5minstr, bin5minquant];

 
% Para 3 min:

bin3minquant = cell2mat(Data83(:, 2:end)); bin3minquant = array2table(bin3minquant, 'VariableNames', ["BinDist 3min Time (s)", "BinMean 3min Speed (cm/s)", "BinDistance 3min (cm)", "BinFreezEpisodes 3min (n)", "BinTotalDurEpi 3min (s)", "BinMeanDurEpi 3min (s)", "BinFreezeEpiNorm 3min (s/min)", "BinTotalDurepiNorm 3min (s/min)"]);
bin3minstr = repmat(catnames.', [length(videopaths) 1]); bin3minstr = array2table(bin3minstr, 'VariableNames', "Bin Number");
mdbdata3 = [bin3minstr, bin3minquant];


 
% % De variables independientes del freezing:
% 
% fig1 = figure(1)
% tiledlayout(3, 1)
% 
% nexttile
%         boxchart(mdbdata.bin_number, mdbdata.time_in_bin, 'BoxFaceColor', [0.4940 0.1840 0.5560], ...
%         'MarkerStyle', '+', 'MarkerSize', 4, 'MarkerColor', 'black')
%         hold on
%         plot(mean(mdbdata.time_in_bin, 1), 'LineStyle', '--', 'LineWidth', 0.5, 'Color', 'black', 'Marker', '.')
%         title('Min Distance to Border vs Time')
%         hold off
% 
% nexttile
%         boxchart(mdbdata.bin_number, mdbdata.mean_speed, 'BoxFaceColor', [0.4940 0.1840 0.5560], ...
%         'MarkerStyle', '+', 'MarkerSize', 4, 'MarkerColor', 'black')
%         hold on
%         plot(mean(mdbdata.mean_speed, 1), 'LineStyle', '--', 'LineWidth', 0.5, 'Color', 'black', 'Marker', '.')
%         title('Min Distance to Border vs Mean Speed')
%         hold off
% 
% nexttile
%         boxchart(mdbdata.bin_number, mdbdata.distance_in_bin, 'BoxFaceColor', [0.4940 0.1840 0.5560], ...
%         'MarkerStyle', '+', 'MarkerSize', 4, 'MarkerColor', 'black')
%         hold on
%         plot(mean(mdbdata.distance_in_bin, 1), 'LineStyle', '--', 'LineWidth', 0.5, 'Color', 'black', 'Marker', '.')
%         title('Min Distance to Border vs Distance Traveled')
%         hold off
% 
% % Variables dependientes del freezing sin normalizar:
% 
% fig2 = figure(2)
% tiledlayout(3, 1)
% 
% nexttile
%         boxchart(mdbdata.bin_number, mdbdata.num_of_epis, 'BoxFaceColor', [0.4940 0.1840 0.5560], ...
%         'MarkerStyle', '+', 'MarkerSize', 4, 'MarkerColor', 'black')
%         hold on
%         plot(mean(mdbdata.num_of_epis, 1), 'LineStyle', '--', 'LineWidth', 0.5, 'Color', 'black', 'Marker', '.')
%         title('Min Distance to Border vs Number of Episodes')
%         hold off
% 
% nexttile
%         boxchart(mdbdata.bin_number, mdbdata.total_dur, 'BoxFaceColor', [0.4940 0.1840 0.5560], ...
%         'MarkerStyle', '+', 'MarkerSize', 4, 'MarkerColor', 'black')
%         hold on
%         plot(mean(mdbdata.total_dur, 1), 'LineStyle', '--', 'LineWidth', 0.5, 'Color', 'black', 'Marker', '.')
%         title('Min Distance to Border vs Total Duration of Episodes')
%         hold off
% 
% nexttile
%         boxchart(mdbdata.bin_number, mdbdata.mean_dur, 'BoxFaceColor', [0.4940 0.1840 0.5560], ...
%         'MarkerStyle', '+', 'MarkerSize', 4, 'MarkerColor', 'black')
%         hold on
%         plot(mean(mdbdata.total_dur, 1), 'LineStyle', '--', 'LineWidth', 0.5, 'Color', 'black', 'Marker', '.')
%         title('Min Distance to Border vs Mean Duration of Freezing')
%         hold off
% 
% 
% 
% % Variables dependientes del freezing normalizadas:
% 
% 
% fig3 = figure(3)
% tiledlayout(2, 1)
% 
% nexttile
%         boxchart(mdbdata.bin_number, mdbdata.freez_epi_norm, 'BoxFaceColor', [0.4940 0.1840 0.5560], ...
%         'MarkerStyle', '+', 'MarkerSize', 4, 'MarkerColor', 'black')
%         hold on
%         plot(mean(mdbdata.freez_epi_norm, 1), 'LineStyle', '--', 'LineWidth', 0.5, 'Color', 'black', 'Marker', '.')
%         title('Min Distance to Border vs Number of Episodes (%)')
%         hold off
% nexttile
%          boxchart(mdbdata.bin_number, mdbdata.dur_epi_norm, 'BoxFaceColor', [0.4940 0.1840 0.5560], ...
%         'MarkerStyle', '+', 'MarkerSize', 4, 'MarkerColor', 'black')
%         hold on
%         plot(mean(mdbdata.dur_epi_norm, 1), 'LineStyle', '--', 'LineWidth', 0.5, 'Color', 'black', 'Marker', '.')
%         title('Min Distance to Border vs Total Duration of Episodes (%)')
%         hold off
% 
% 
% % Guardar datos en formato xlsx y figuras en formato svg
% 
 xlsx8 = strcat(metricspath, '/MDBxContVars.xlsx');
 xlsx833 = strcat(metricspath, '/MDBxContVars3min.xlsx');
% figpath1 = strcat(metricspath, '/MDBxIndofFreez.svg');
% figpath2 = strcat(metricspath, '/MDBxDepofFreez.svg');
% figpath3 = strcat(metricspath, '/MDBxDepofFreezNorm.svg');
% figpath4 = strcat(metricspath, '/MDBxIndofFreez.jpeg');
% figpath5 = strcat(metricspath, '/MDBxDepofFreez.jpeg');
% figpath6 = strcat(metricspath, '/MDBxDepofFreezNorm.jpeg');
% 





 writetable(mdbdata, xlsx8, 'WriteVariableNames', true)
writetable(mdbdata3, xlsx833, 'WriteVariableNames', true)
% saveas(fig1, figpath1)
% saveas(fig2, figpath2)
% saveas(fig3, figpath3)
% saveas(fig1, figpath4)
% saveas(fig2, figpath5)
% saveas(fig3, figpath6)

end % Fin de la función 
