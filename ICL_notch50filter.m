
function [] = ICL_notch50filter(folderpath, varargin)

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
for i = 1:length(listOfFolderNames)
%     if contains(listOfFolderNames(i), 'Data') && contains(listOfFolderNames(i), 'continuous')
    if contains(listOfFolderNames(i), 'recording') && ~contains(listOfFolderNames(i), 'events') && ~contains(listOfFolderNames(i), 'continuous')
        lfp_paths = [lfp_paths listOfFolderNames(i)];
    else
        continue
    end
end


%% 2. APPLY FILTER:

for ii = 1:length(lfp_paths)
    
    % 2.1. Load file:

    rawname = char(strcat(lfp_paths(ii),'\LFPraw.mat'));
    resultsname = char(strcat(lfp_paths(ii), '\LFPdsdata.mat'));
    load(rawname,'results_rawLFP');
    load(resultsname, 'results');

        message = sprintf(['Processing file %i out of %i' ...
        ' ' ...
        ], ii, length(lfp_paths));
        disp(message)

    
  
    % 2.2. Load variables:
    
    LFPdata = results_rawLFP.LFP.LFPdata;
    fs = results_rawLFP.LFP.fs;


    % CÓDIGO DE ALBERTO: Quitar ruido 50 Hz (notch filter)
    wo = 50/(fs/2);   factor  = 10;
    bw= wo/ factor;
    [b, a] = iirnotch(wo, bw);
    data_filt = filtfilt(b,a,LFPdata);
    LFPdata = data_filt;

    % 2.3. Downsample original data with notch filter in order to save in
    % dsdata for further analysis:

    newfs = dsfs;
    fsdown = fs./newfs;
    dsdata = downsample(LFPdata, fsdown);


    % 2.4. Save variables to LFPdsdata file:

    results.LFPds.dsdata = dsdata;
    results.LFPds.fs = newfs;
    
    % 2.5. Save file:

    save(resultsname, 'results')

end


end