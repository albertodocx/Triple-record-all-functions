
%%% Función para pasar datos de BHD a CSV: (Para lo de Nuria)


%%% BHDextract.mat extracts relevant variables for posterior analysis:

% Event and timing of event (freezing, licking, shock, etc)
% Location
% Velocity 

%clear all; clc;

function [] = BHDtoXLSXbatch2()

% Prueba para coger muchos archivos: 

path = uigetdir();

%%%% ESTO ES DE IMAGEANALYST MATHWORKS:

% Get list of all subfolders.
allSubFolders = genpath(path);

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

%%%%% hasta aquí lo de imageanalyst mathworks

videopaths = [];
for i = 1:length(listOfFolderNames)
    if contains(listOfFolderNames(i), 'analyzed')
        videopaths = [videopaths listOfFolderNames(i)];
    else
        continue
    end
end


for i = 1:length(videopaths)

BHDpath = string(videopaths(i));

% Get path for BHD *.mat files. 

Behpath = strcat(BHDpath, '/Behavior.mat');
Metpath = strcat(BHDpath, '/Metrics.mat');
Fspath = strcat(BHDpath, '/Params.mat');

% Load files*.mat to function workspace 

load(Behpath);
load(Metpath);
load(Fspath);

% Save variables to struct results. 

BHD = struct();
BHD.Event.Time = Behavior.Freezing.Bouts;
BHD.Event.Type = 'Freezing';
BHD.Location = round(Metrics.Location.'); % La transpuesta porque está el eje x en una fila y el eje y en otra fila
BHD.Fs = Params.Video.frameRate;
BHD.MidbackSpeed = round(Metrics.Velocity.MidBack.', 3); % Esta variable Midback y Movement.Data dan el mismo vector, la diferencia es que Midback está corregido para tener 4500 muestras y no n-1 que resulta de calcular diff. 
BHD.path = BHDpath;
BHD.DistanceTraveled = Metrics.Movement.DistanceTraveled.';


% Sacar vector lógico de freezing: 

BHD.freezlogical = Behavior.Freezing.Vector.';


% Crear vector de tiempo: 

Time = round(((0:length(BHD.Location)-1) ./ BHD.Fs), 5);
Time = Time.';
BHD.Time = Time; % Guardarlo en struct BHD para guardar .mat con mismos datos que CSV


% Calcular velocidad:


vel = diff(BHD.Location) ./ BHD.Fs; %cada fila de la posicion es un punto (tantas columnas como ejes)

vel = [vel(1,:); vel]; % Repite el primer valor puesto que diff saca length-1

vel = sqrt(vel(:,1).^2 + vel(:,2).^2); % Puesto que Nuria tiene dos ejes, hay que sacar la velocidad global (el módulo de dos vectores)

BHD.Velocity = vel;

% Complete Data: 

Data = [Time, BHD.Location, vel, BHD.MidbackSpeed, BHD.freezlogical];

colnames = {'Time', 'X axis', 'Y axis', 'Global Velocity', 'BHD Velocity', 'Freezing (log)'};

% Convertir a tabla para poder guardar en csv..? 

Data = array2table(Data, 'VariableNames', colnames);

% Guardar en csv y en .mat

savepath = strcat(BHDpath, '/BehData.xlsx'); % Mirar de cambiar a CSV para que se pueda importar a otros programas 
writetable(Data, savepath);
%csvwrite(savepath, Data);

savemat = strcat(BHDpath, '/BehData.mat');
save(savemat, 'BHD');

end

end



