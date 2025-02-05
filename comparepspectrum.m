
function [] = comparepspectrum()

%% 0. Automatic loading of all files in a folder: 

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

powfqdata = zeros(1999, length(lfp_paths) + 1);
orderinfo = array2table(strings(length(lfp_paths), 3), 'VariableNames', ["mouseid", "genotype", "path"]);
orderinfo.mouseid = categorical(orderinfo.mouseid);
orderinfo.genotype = categorical(orderinfo.genotype);
orderinfo.path = categorical(orderinfo.path);

for nn = 1:length(lfp_paths)

    %% 1. Load saved file and other variables:

    filename = char(strcat(lfp_paths(nn), '\LFPdsdata.mat'));
    message = sprintf(['Processing file %i out of %i' ...
        ' ' ...
        ], nn, length(lfp_paths));
    disp(message)
    load(filename, 'results')
%% 2. 
    psProfile = results.psProfile.pS;
    pSfreqs = results.psProfile.pSfreqs; % x axis // Fq(Hz)

    meanpower = pow2db(mean(mean(psProfile, 2), 3));
    meanpower(pSfreqs > 45 & pSfreqs < 55) = nan;  
    meanpower = fillmissing(meanpower, 'linear');

    figure(nn)
    plot(pSfreqs, meanpower)
    xlabel('Frequency (Hz)')
    ylabel('Power (dB)')


    %     meanpowerint = interp1(pSfreqs, meanpower, pSfreqs(pSfreqs > 45 & pSfreqs < 55));
    %     meanpower(pSfreqs > 45 & pSfreqs < 55) = meanpowerint;

    powfqdata(:, nn+1) = meanpower;

    orderinfo.mouseid(nn, 1) = results.id.mouseid;
    orderinfo.genotype(nn, 2) = results.id.expcond;
    orderinfo.path(nn, 3) = results.path;
end

powfqdata(:, 1) = pSfreqs;


       savename = uigetdir();
       dataname = strcat(savename, '\compare_pssummary.xlsx');
       dataname2 = strcat(savename, '\looporderinfo.xlsx');
       
       writematrix(powfqdata, dataname);
       writetable(orderinfo, dataname2, 'WriteVariableNames', true)

end