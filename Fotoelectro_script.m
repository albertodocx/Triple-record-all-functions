%% Script para analizar fotometría,electrofisiología y conducta
% Las funciones en la que me he ido basando son anteriores, de Rut

%% 1. CARGAR CARPETA Y SUBCARPETAS CON LOS ARCHIVOS BRUTOS
% Carpeta principal
mainFolder = uigetdir;

% Obtener lista de subcarpetast
subfolders = dir(mainFolder);
subfolders = subfolders([subfolders(:).isdir] & ~ismember({subfolders(:).name}, {'.', '..'}));

%% Reordenar


%% 2. ABRIR LOS DATOS BRUTOS EN FORMATO TDTdata
% Iterar a través de las subcarpetas
for i = 1:length(subfolders)
    disp(['vamos por i = ', num2str(i)]);
    currentSubfolder = fullfile(mainFolder, subfolders(i).name);
    
    % Ejecutar el script en la subcarpeta actual
    mypath = currentSubfolder;
    filename = 'fpdata.mat';
    TDTdata = fp2mat(mypath, filename);
    resultsname = strcat(currentSubfolder,'/fpdata.mat');
    save(resultsname, 'TDTdata');
end

%% 3. PREPROCESAMIENTO DE SEÑALES Y AJUSTE A CONDUCTA
for i = 1:length(subfolders)
    disp(['vamos por i = ', num2str(i)]);
    currentSubfolder = fullfile(mainFolder, subfolders(i).name);
    mypath = currentSubfolder;
    mypathTDT = strcat(currentSubfolder,'/fpdata.mat');
    load (mypathTDT,"TDTdata")

    TDTdata.currentsignals = conducta (TDTdata);
    resultsname = strcat(currentSubfolder,'/fpdata.mat');
    save(resultsname, 'TDTdata');
    close all
end

%% 4. ANALIZAR FOTOMETRÍA 
for i = 1:length(subfolders)
    disp(['vamos por i = ', num2str(i)]);
    currentSubfolder = fullfile(mainFolder, subfolders(i).name);
    mypath = currentSubfolder;
    mypathTDT = strcat(currentSubfolder,'/fpdata.mat');
    load (mypathTDT,"TDTdata") 

    % Datos de fotometría
    results = fp_preproc_Alberto(TDTdata,'SampleBuffer',[0, 0]);
    %cortarpresentacion(results, mypath);%cortar gráfica
    close all
end

%% 5. ANALIZAR ELECTROFISIOLOGÍA
for i = 1:length(subfolders)
    disp(['vamos por i = ', num2str(i)]);
    currentSubfolder = fullfile(mainFolder, subfolders(i).name);
    mypath = currentSubfolder;
    mypathTDT = strcat(currentSubfolder,'/fpdata.mat');
    mypathRs = strcat(currentSubfolder,'/Results.mat');
    load (mypathTDT,"TDTdata")
    load (mypathRs,"results")

 % Datos de electrofisiología

    datosFiltrados = notchfilter(TDTdata);
    results.lfp.datosfiltrados = datosFiltrados;
    potenciasLFP(TDTdata, results, mypath);
    
    close all
    
    pyr = results.lfp.datosfiltrados.';
    pyrfs = TDTdata.streams.Wav1.fs;
    pyrdata.dataLFP = pyr;
    pyrdata.fsLFP = pyrfs;

    psprofilebueno = LFP_compute_psProfile(pyrdata);
    results.lfp.psprofilebueno = psprofilebueno; %revisar para unir esto con lo anterior
    close all

    espectrograma(TDTdata,results,mypath);
    close all

    results.lfp.thetaanalisis = theta_eventsAlb(results,mypath);
    close all

    resultsname = strcat(currentSubfolder,'/Results.mat');
    save(resultsname, 'results');

    results.lfp.lfp_zscore = lfpzscore(results,mypath);
    
    resultsname = strcat(currentSubfolder,'/Results.mat');
    save(resultsname, 'results');

end

%% 6. ANALIZAR FOTOMETRÍA JUNTO CON ELECTROFISIOLOGÍA

for i = 1:length(subfolders)
    disp(['vamos por i = ', num2str(i)]);
    currentSubfolder = fullfile(mainFolder, subfolders(i).name);
    mypath = currentSubfolder;
    mypathTDT = strcat(currentSubfolder,'/fpdata.mat');
    mypathRs = strcat(currentSubfolder,'/Results.mat');
    load (mypathTDT,"TDTdata")
    load (mypathRs,"results")

 % BEH-PETH análisis LFP y fotometría en conjunto
 
    DFF = results.FP.Signals.DFFModZscore;
    fsDFF = TDTdata.streams.x465N.fs;
    eventos = floor(results.lfp.thetaanalisis.bouts.muestras);%revisar (2)
    fsev = 2;%revisar (2) 
    maxlength = floor(length(DFF)/fsDFF);


    results.FP.analysis = BEH_PETHonsetv2(DFF, fsDFF, eventos, fsev, 'dffmovement', mypath,'maxlength',maxlength,'AUCint',[-7 0; 0 7], 'Pre', 7, 'Post', 7, 'bin', 0.2);
    %results.ripplesresults = getripplesAlb(TDTdata,results,mypath);

    resultsname = strcat(currentSubfolder,'/Results.mat');
    save(resultsname, 'results');
    
end

%% 7. ANALIZAR CONDUCTA
for i = 1:length(subfolders)
    disp(['vamos por i = ', num2str(i)]);
    currentSubfolder = fullfile(mainFolder, subfolders(i).name);
    mypath = currentSubfolder;
    mypathTDT = strcat(currentSubfolder,'/fpdata.mat');
    mypathRs = strcat(currentSubfolder,'/Results.mat');
    load (mypathTDT,"TDTdata")
    load (mypathRs,"results")
    
    % results.DLC = openDLC(mypath);
    % close all
    % results.DLC.bouts = dlcmovement(results.DLC.data.manos,mypath);
    % close all
    % resultsname = strcat(currentSubfolder,'/Results.mat');
    % save(resultsname, 'results');
    % 
    % % BEH-PETH análisis conducta y fotometría en conjunto
    % DFF = results.FP.Signals.DFFModZscore;
    % fsDFF = TDTdata.streams.x465N.fs;
    % eventos = floor((results.DLC.bouts.segcut)*15);
    % fsev = 15;
    % maxlength = floor(length(DFF)/fsDFF);
    % path = strcat(mypath,'/conducta');
    % 
    % if isnan(eventos)
    %     disp ('otra vez será')
    % else
    % 
    % results.DLC.conducta = BEH_PETHonsetv2(DFF, fsDFF, eventos, fsev, 'dffmovement', path,'maxlength',maxlength,'AUCint',[-7 0; 0 7], 'Pre', 7, 'Post', 7, 'bin', 0.2);
    % 
    % resultsname = strcat(currentSubfolder,'/Results.mat');
    % save(resultsname, 'results');
    % end

    % BEH-PETH análisis conducta y theta en conjunto
    DFF = results.lfp.lfp_zscore.data.thetazscore;
    fsDFF = 2;
    eventos = floor((results.DLC.bouts.segcut)*15);
    fsev = 15;
    maxlength = floor(length(DFF)/fsDFF);
    path = strcat(mypath,'/conductath');

    if isnan(eventos)
        disp ('otra vez será')
    else
    results.DLC.conductath = BEH_PETHonsetv2(DFF, fsDFF, eventos, fsev, 'dffmovement', path,'maxlength',maxlength,'AUCint',[-7 0; 0 7], 'Pre', 7, 'Post', 7, 'bin', 0.2);
    resultsname = strcat(currentSubfolder,'/Results.mat');
    save(resultsname, 'results');
    end

end

%% Análisis completo final:

[~, nombre_raton, ~] = fileparts(subfolders(1).name); %Cuidado, si hay más de un ratón en la carpeta solo cogerá 1 de ellos para esta variable
idx_medio_guion = regexp(nombre_raton, '(?<=_)[A-Za-z]\d+(?=_)', 'match', 'once');
nombre_raton = idx_medio_guion;


if ~exist(fullfile(mainFolder, ['Análisis final ',nombre_raton]), 'dir')
    mkdir(fullfile(mainFolder, ['Análisis final ',nombre_raton]));
end

% Convertir nombre_raton a una matriz de caracteres de una sola fila
nombre_raton = char(nombre_raton);

finalfolder = fullfile(mainFolder, ['Análisis final ', nombre_raton]);


%savefinal (mainFolder,subfolders,finalfolder);

% saveanimalpower (mainFolder,subfolders,finalfolder);

% saveanalisisFPLFP (mainFolder,subfolders,finalfolder);
% saveanalisisFPCond (mainFolder,subfolders,finalfolder);
% saveanalisisthetaLFP (mainFolder,subfolders,finalfolder);
saveanalisisthetaCond (mainFolder,subfolders,finalfolder);

