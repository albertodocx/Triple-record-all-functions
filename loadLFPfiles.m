 
function [] = loadLFPfiles()

% Load main folder: 

uiwait(msgbox('Select main path (main folder)', 'Instructions', "modal")); 

folderpath = uigetdir(); 


%% 1. BATCH PROCESSING:

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
%     if contains(listOfFolderNames(i), 'Data') & contains(listOfFolderNames(i), 'continuous')
    if contains(listOfFolderNames(i), 'recording') && ~contains(listOfFolderNames(i), 'events') && ~contains(listOfFolderNames(i), 'continuous')
        lfp_paths = [lfp_paths listOfFolderNames(i)];
    else
        continue
    end
end

%% 2. LOOP THROUGH ALL FILES:

for n = 1:length(lfp_paths)

%     filename = char(strcat(lfp_paths(n), '\continuous.dat'));
    filename = char(strcat(lfp_paths(n), '\structure.oebin'));
    message = sprintf(['Loading file %i out of %i' ...
        ' ' ...
        ], n, length(lfp_paths));
    disp(message)

    %%  3. Open file:
% 
%     f=fopen(filename,'rb');
%     nChannels = 32;
%     D=fread(f,[nChannels, Inf],'int16');
%     fclose(f);

     data = load_open_ephys_binary(filename, "continuous", 1);

    %%  4. Get only 16 channels (rest not used):

    LFPdata = data.Data(1:16, :);
    bitvolts = data.Header.channels.bit_volts;
    LFPdata = LFPdata * bitvolts;

    fs = 20000; % Hz
    lengthinsec = size(LFPdata, 2)/fs;
    lengthinmin = lengthinsec/60; 
    LFPdata = LFPdata.';

    %% 5.  Save LFPdata in raw version and another file for dsdata:

    results_rawLFP.LFP.LFPdata = LFPdata;
    results_rawLFP.LFP.fs = fs;
    results_rawLFP.LFP.durinsec = lengthinsec;
    results_rawLFP.LFP.durinmin = lengthinmin;
    results_rawLFP.path = filename;

    results.params.fs = fs;
    results.params.durinsec = lengthinsec;
    results.params.durinmin = lengthinmin;
    results.params.path = filename;


    
    %% 6. Downsample data: 

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

%     for i = 1:ch
%         dsdata(:, i) = downsample(LFPdata(:, i), fsdown);
%     end

    results.LFPds.dsdata = dsdata;
    results.LFPds.fs = newfs;

    rawname = char(strcat(lfp_paths(n), '\LFPraw.mat'));
    resultsname = char(strcat(lfp_paths(n), '\LFPdsdata.mat'));
    save(rawname, 'results_rawLFP', '-v7.3') % Estos datos tardan mucho en cargar y por eso guardo otro archivo que sea más manejable posteriormente
 %   save(resultsname, 'results', '-v7.3')
    save(resultsname, 'results')

end % Bucle de todos los archivos

end % Fin de la función
