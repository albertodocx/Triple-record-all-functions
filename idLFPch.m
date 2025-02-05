
function [] = idLFPch()

%% 1. Open and save folder:

% uiopen()
% 
% savepath = results.params.path;
% savepath = erase(savepath, "\continuous.dat");
% folderpath = savepath;

%%%% ESTO ES DE IMAGEANALYST MATHWORKS:
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
    if contains(listOfFolderNames(i), 'recording') && ~contains(listOfFolderNames(i), 'events') && ~contains(listOfFolderNames(i), 'continuous')
        lfp_paths = [lfp_paths listOfFolderNames(i)];
    else
        continue
    end
end


for n = 1:length(lfp_paths)

    filename = char(strcat(lfp_paths(n), '\LFPdsdata.mat'));
    message = sprintf(['Processing file %i out of %i' ...
        ' ' ...
        ], n, length(lfp_paths));
    disp(message)
    load(filename)

% 2. Check which channels to keep with power spectrum:

% 2.1. Power spectrum:

%ch = size(results.LFP.LFPdata, 1); Siempre serán 16 canales:

ch = 16;

LFPdata = results.LFPds.dsdata;
fs = results.LFPds.fs;

 psProfile = LCN_compute_psProfile(LFPdata, fs, ch);
%   psProfile = powerSpectrumProfileAlb(LFPdata, fs, 1);

% 2.2. Select channels of interest (theta + ripples):

 prompt = {'Ripples Channel', 'Theta Channel'};

 dlgtitle = 'Select channels of interest';
 dims = [1 35];
 definput = {'4', '16'};
 
 SB = inputdlg(prompt, dlgtitle, dims, definput);
 
 ripplech = str2num(SB{1});
 thetach = str2num(SB{2});
 
 ripple_chname = sprintf('ch%i', ripplech);
 theta_chname = sprintf('ch%i', thetach);

% 2.3. Introduce labels to identify subject, condition and experiment:

prompt = {'filename', 'Mice ID:','Experimental Condition', 'Pre vs Post'};
dlgtitle = 'Enter labels for file';
dims = [4 40];
definput = {char(lfp_paths(n)),'C101 44B2','VCD', 'pre'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

% 3.  Save all variables:

% 3.1. Save figs:

figname = strcat(lfp_paths(n), '\psperch.jpg');
fig = figure(n);
saveas(fig, char(figname));

% 3.2. Save vars:

results.folderpath = folderpath;
results.params.ripplech = ripplech;
results.params.thetach = thetach;
% results.params.ripplech = ripch;
% results.params.thetach = thch;
results.psProfile = psProfile;
results.id.filename = answer{1};
results.id.mouseid = answer{2};
results.id.expcond = answer{3};
results.id.time = answer{4};


resultsname = char(strcat(lfp_paths(n), '\LFPdsdata.mat'));
results.path = resultsname;
save(resultsname, 'results')
%save(resultsname, 'results', '-v7.3')

end % del bucle
end % de la función

