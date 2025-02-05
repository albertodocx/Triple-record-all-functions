
%%% Shockmeasuresloc.m extracts the following measures for experiments
%%% where the shock is given depending on animal's location (center/periphery)

%%% De reactancia: 

% 1. Position at shock time (tiempo total - 30 segs) >> pendiente de que
% Nuria confirme
% 2. Mean speed after shock ([28 y 26 seg antes de que termine el registro])
% 20.04.2023: visto con Nuria, mean speed after shock durante el shock de
% ttotal -30 a -28 seg (y no lo que tenía puesto de 28 y 26). 

%%% En relación con freezing:

% 3. Dur total del freezing después del shock
% 4. Dur media del freezing después del shock
% 5. Num episodios después  del shock
% 6. Dur cada episodio después del shock 


function [] = shockmeasuresloc()

% 1. Select folder where files will be saved:

uiwait(msgbox('Select path (a folder/directorio) where you want to save your complete data', 'Instructions', "modal")); 

metricspath = uigetdir(); % Se utiliza esta ruta para generar los tres *.mat. 

% 2. Select main folder to analyze shocks

    message = 'Please select main folder to analyze';
    uiwait(msgbox(message, 'Instructions', 'modal'));

    mainfolder = uigetdir()

% Get list of all subfolders.
allSubFolders = genpath(mainfolder);

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


shockpaths = [];
for i = 1:length(listOfFolderNames)
    if (contains(listOfFolderNames(i), 'Day02') || contains(listOfFolderNames(i), 'Day2')) & (contains(listOfFolderNames(i), 'analyzed'))
        shockpaths = [shockpaths listOfFolderNames(i)];
    else
        continue
    end
end


% Create datasheet:

Data = {};
filename = strcat(metricspath, '/CenterPeriShockMeasures.mat');
measurenames = { 'Test Duration (s)', 'Shock Event (s)' ,'Shock Location in X','Shock Location in Y', ...
                 'Mean Vel during Shock (cm/s)','Min Dist at Shock 0s (cm)', 'Min Dist after Shock 2s (cm)', 'path'};


for i = 1:length(shockpaths)
    
% Select animal
    subpath = strcat(string(shockpaths(i)), '\BehData.mat');
    load(subpath)

% Extract variables


    %distBorder = BHD.distToBorder;
    Location = BHD.Location;
    Fs = BHD.Fs;

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


% Create and calculate variables:

    % X,Y, coordinates at the time of the shocks:
    
    testdur = length(Location)/Fs;
    shocktime = length(Location) - round(30*Fs);
    shockevent = shocktime/Fs;

    loc_shocktime = BHD.Location(shocktime,:);
    
    % Distance to closest border at the time of the shocks:

    distborder_shocktime = minDist(shocktime);
  
   % Distance to closest border immediately after the shocks:
    
    distborder_aftershock = minDist(shocktime + round(2*Fs));

    % Velocidad media: 

    meanvel_duringshock = mean(BHD.MidbackSpeed(shocktime:shocktime + round(2*Fs)));
    

% Save datA:

Data(i, :) = {testdur, shockevent, loc_shocktime(1), loc_shocktime(2),... 
        meanvel_duringshock, distborder_shocktime, distborder_aftershock, subpath}; 

save(filename, 'Data')


end % fin del bucle

    Data = array2table(Data, 'VariableNames', measurenames);
    xlsx = strcat(metricspath, '/CenterPeriShockMeasures.xlsx');
    writetable(Data, xlsx, 'WriteVariableNames', true); 



end % de la función