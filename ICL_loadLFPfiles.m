 
function [] = ICL_loadLFPfiles(folderpath, varargin)

% INPUTS:
%  folderpath: main folder where data is stored in subfolders
%  dsfs: new frequency sample for downsampling. Default is 2000 Hz. 
% NO OUTPUTS (files are saved in each path)

p = inputParser;
addParameter(p,'dsfs', 2000 ,@isnumeric);
parse(p,varargin{:});
dsfs = p.Results.dsfs;


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
for ii = 1:length(listOfFolderNames)
%     if contains(listOfFolderNames(i), 'Data') & contains(listOfFolderNames(i), 'continuous')
    if contains(listOfFolderNames(ii), 'recording') && ~contains(listOfFolderNames(ii), 'events') && ~contains(listOfFolderNames(ii), 'continuous')
        lfp_paths = [lfp_paths listOfFolderNames(ii)];
    else
        continue
    end
end

%% 2. LOOP THROUGH ALL FILES:

for nn = 1:length(lfp_paths)

    filename = char(strcat(lfp_paths(nn), '\structure.oebin'));
    message = sprintf(['Loading file %i out of %i' ...
        ' ' ...
        ], nn, length(lfp_paths));
    disp(message)

    %%  3. Open file:

     data = load_open_ephys_binary(filename, "continuous", 1);

    %%  4. Get only 16 channels (rest not used):

    LFPdata = data.Data(1:16, :);
    bitvolts = data.Header.channels.bit_volts;
    LFPdata = LFPdata * bitvolts;

    fs = data.Header.sample_rate;
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
    
    newfs = dsfs;
    fsdown = fs./newfs;

    disp('downsampling...'); dsdata = downsample(LFPdata,fsdown);

    results.LFPds.dsdata = dsdata;
    results.LFPds.fs = newfs;

    %% 7. Save data: 

    rawname = char(strcat(lfp_paths(nn), '\LFPraw.mat'));
    resultsname = char(strcat(lfp_paths(nn), '\LFPdsdata.mat'));
    save(rawname, 'results_rawLFP', '-v7.3') 
    save(resultsname, 'results')

end % Files

end % Function
