
%%% Shockmeasures.m extracts the following measures for experiments with 3
%%% standardized shocks at 180 s, 210 s and 240 s:

% 1. Position at 180s - xy
% 2. Position at 210s - xy
% 3. Position at 240s - xy
% 4. Mean speed after first shock (182-184)
% 5. Mean speed after second shock (212-214)
% 6. Mean speed after third shock (240-242)
% 7. Mean of speed means (4,5,6). 

% 20.04.2023: cambiar 4, 5 y 6 a 180-182/210-212/240-242


function [] = shockmeasures()

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
    filename = strcat(metricspath, '/ShockMeasures.mat');
    measurenames = { '180 X', '180 Y', ...
                 '210 X','210 Y',...
                 '240 X', '240 Y',...
                 'Mean Vel 180 2 s (cm/s)', 'Mean Vel 210 2 s (cm/s)', 'Mean Vel 240 2 s (cm/s)', 'Total Mean Vel (cm/s)',...
                 'Min Dist 180 (cm)', 'Min Dist 210 (cm)', 'Min Dist 240 (cm)', ...
                 'Min Dist 182 (cm)', 'Min Dist 212 (cm)', 'Min Dist 242 (cm)', 'path' };


% Load each file:


for i = 1:length(shockpaths)

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

    loc180 = BHD.Location(round(180*Fs),:);
    loc210 = BHD.Location(round(210*Fs),:);
    loc240 = BHD.Location(round(240*Fs),:);
    
    % Distance to closest border at the time of the shocks:

    distborder180 = minDist(round(180*Fs)); 
    distborder210 = minDist(round(210*Fs));
    distborder240 = minDist(round(240*Fs));
  
   % Distance to closest border immediately after the shocks:
    
    distborder182 = minDist(round(182*Fs));
    distborder212 = minDist(round(212*Fs));
    distborder242 = minDist(round(242*Fs));


    % Velocidad media: 

    meanvel180 = mean(BHD.MidbackSpeed(round(180*Fs):round(182*Fs)));
    meanvel210 = mean(BHD.MidbackSpeed(round(210*Fs):round(212*Fs)));
    meanvel240 = mean(BHD.MidbackSpeed(round(240*Fs):round(242*Fs)));
    
    totalvelduringshock = mean([meanvel180 meanvel210 meanvel240]);


% Save datA:

    Data(i, :) = {loc180(1, 1), loc180(1, 2),...
                loc210(1, 1), loc210(1, 2),...
                loc240(1, 1), loc240(1, 2),...
                meanvel180, meanvel210, meanvel240, totalvelduringshock,...
                distborder180, distborder210, distborder240,...
                distborder182, distborder212, distborder242, subpath};


    save(filename, 'Data')
    

end % fin del bucle

    Data = array2table(Data, 'VariableNames', measurenames);
    xlsx = strcat(metricspath, '\ShockMeasures.xlsx');
    writetable(Data, xlsx, 'WriteVariableNames', true); 



end % de la función