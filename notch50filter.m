
function [] = notch50filter()

%% 1. BATCH PROCESSING:

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


%% 2. APPLY FILTER:

for ii = 1:length(lfp_paths)
    
    % 2.1. Load file:

    rawname = char(strcat(lfp_paths(ii),'\LFPraw.mat'));
    resultsname = char(strcat(lfp_paths(ii), '\LFPdsdata.mat'));
    load(rawname);
    load(resultsname);

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

    ch = size(LFPdata, 2);
   
    lendata = zeros(1, size(LFPdata, 1));
    fsdown = fs./2000;
    newfs = 2000;
    lendata = length(downsample(lendata, fsdown));
    dsdata = zeros(lendata, ch);
    dsdata = downsample(LFPdata, fsdown);


    % 2.4. Save variables to LFPdsdata file:

    results.LFPds.dsdata = dsdata;
    results.LFPds.fs = newfs;
    
    % 2.5. Save file:

    save(resultsname, 'results')

end


end