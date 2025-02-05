
function alldata = BEH_PHOT_binvar(contvar, DFF,behfs,fs, threshold, binsize, savepath)

% INPUTS:
%  DFF: continuous variable with DFF photometry data. Use dFF Z-Score if
% available. 
%  fs: fiber photometry frequency sampling 
%  contvar: continuous variable to correlate with DFF
%  behfs: behavioral frequency sampling (video, etc) 
%  threshold: a vector with two values; min and max of bins. 
%  binsize: how many bins to use

% OUTPUTS:
%  bins: vector with the values of the bins

%% 1. Create timestamps for each variable: 

ts = (0:length(DFF)-1)./fs;
behts = (0:length(contvar)-1)./behfs;

% 1.1. FF signal has usually higher fs: downsample to ease the calc 

dsfactor = round(fs / behfs);
DFFds = downsample(DFF, dsfactor);
tsds = (0:length(DFFds)-1)./behfs;

% 1.2. Smooth continuous signal (behfs is the smoothing factor):

contvar = smooth(contvar, behfs, 'moving');

%% 2. Create bin vector: 

bins = linspace(threshold(1), threshold(2), binsize+1); % 11 values for 10 bins

% 2.1. Preallocate each bin as a logical matrix: 

binidx = false(binsize, length(DFFds));

% Correct for different lengths:
if length(DFFds) < length(contvar)
    contvar = contvar(1:length(DFFds));
elseif length(DFFds) > length(contvar)
    DFFds = DFFds(1:length(contvar));
end

for ii = 1:binsize
    binidx(ii, :) = contvar>= bins(ii) & contvar < bins(ii+1); 
    binbouts{ii}= findonoff(binidx(ii, :));
end


%% 3. Find onset and offset of each true event within each bin: 

binnedDatads = cell(4, binsize); % 4 measures (signals, AUC, dens, Amplitude - max-)
binnedData = binnedDatads;

meanDFFds = nan(binsize, 1); meanDFF = nan(binsize, 1);
meanAUCds = nan(binsize, 1); meanAUC = nan(binsize, 1);
meanDENSds = nan(binsize, 1); meanDENS =nan(binsize, 1);
meanMAXds = nan(binsize, 1); meanMAX = nan(binsize, 1);


for ii = 1:binsize
    tempbouts = binbouts{1, ii};
    if isempty(tempbouts)
        continue
    else
    tsbouts = [tsds(tempbouts(:, 1)).', tsds(tempbouts(:, 2)).'];
 % Con tsbouts ya se podrían identificar los eventos en la fotometría sin
 % downsamplear
     bDFF = []; bAUC = []; bDens = [];
     bDFFds = []; bAUCds = []; bDensds = [];
     maxDFFds = []; maxDFF = [];

    for jj = 1:size(tempbouts, 1)
        tempDFFdsidx = tsds >= tsbouts(jj, 1) & tsds < tsbouts(jj, 2);
        tempDFFidx = ts >= tsbouts(jj, 1) & ts < tsbouts(jj, 2);     
        tempDFFds = DFFds(tempDFFdsidx);
        tempDFF =  DFF(tempDFFidx);
         if length(tempDFFds) > behfs/2
            tempAUCds = trapz(tsds(tempDFFdsidx), tempDFFds);
            tempDENSds = tempAUCds./length(tempAUCds)./behfs;
         else
            tempAUCds = nan;
            tempDENSds = nan;
         end
         if length(tempDFF) > fs/2
            tempAUC = trapz(ts(tempDFFidx), tempDFF);
            tempDENS = tempAUC./length(tempAUC)./fs;
         else 
             tempAUC = nan;
             tempDENS = nan;
         end
         tempmaxDFFds = max(tempDFFds); tempmaxDFF = max(tempDFF);
         
        % Save variables: 
        bDFFds = [bDFFds, tempDFFds]; bDFF = [bDFF, tempDFF]; % vectores que concatenan vectores 
        bAUCds = [bAUCds, tempAUCds]; bAUC = [bAUC, tempAUC]; % escalar que pasa a vector por cada bout 
        bDensds = [bDensds, tempDENSds]; bDens = [bDens, tempDENS]; % escalar que pasa a vector por cada bout
        maxDFFds = [maxDFFds, tempmaxDFFds]; maxDFF = [maxDFF, tempmaxDFF];
        
    end

    % Means per bin:
    meanDFFds(ii) = mean(bDFFds, 'omitnan'); meanDFF(ii) = mean(bDFF, 'omitnan');
    meanAUCds(ii) = mean(bAUCds, 'omitnan'); meanAUC(ii) = mean(bAUC, 'omitnan');
    meanDENSds(ii) = mean(bDensds, 'omitnan'); meanDENS(ii) = mean(bDens, 'omitnan');
    meanMAXds(ii) = mean(maxDFFds, 'omitnan'); meanMAX(ii) = mean(maxDFF, 'omitnan');

    % Save all the vectors separately per each bin and each variable: 
    binnedDatads{1, ii} = bDFFds; binnedData{1, ii} = bDFF; % DFF Signal
    binnedDatads{2, ii} = bAUCds; binnedData{2, ii} = bAUC; % AUC of signal per each event
    binnedDatads{3, ii} = bDensds; binnedData{3, ii} = bDens; % Density of signal per each event
    binnedDatads{4, ii} = maxDFFds; binnedData{4, ii} = maxDFF;

    end
end

%% 4. Create data tables readable in excel that could be later concatenated with other recordings, files, etc: 

% DFF: 
DFFbineachlen = cumsum(cellfun(@length, binnedData(1, :)));
DFFbinlen = DFFbineachlen(end);
DFFcsvdata = table(strings(DFFbinlen, 1),nan(DFFbinlen, 1), 'VariableNames', ["Bin" "DFF"]);

for ii = 1:binsize
    binvector = repmat(sprintf("Bin%i", ii), length(binnedData{1, ii}), 1);
    if isempty(binvector)
        continue
    else
        if ii == 1
            DFFcsvdata.Bin(1:DFFbineachlen(ii)) = binvector;
            DFFcsvdata.DFF(1:DFFbineachlen(ii)) = binnedData{1, ii}.';
        elseif ii > 1
            DFFcsvdata.Bin(DFFbineachlen(ii-1)+1:DFFbineachlen(ii)) = binvector;
            DFFcsvdata.DFF(DFFbineachlen(ii-1)+1:DFFbineachlen(ii)) = binnedData{1, ii}.';
        end
    end
end


% DFFds: 
DFFdsbineachlen = cumsum(cellfun(@length, binnedDatads(1, :)));
DFFdsbinlen = DFFdsbineachlen(end);
DFFdscsvdata = table(strings(DFFdsbinlen, 1),nan(DFFdsbinlen, 1), 'VariableNames', ["Bin" "DFF"]);

for ii = 1:binsize
    binvector = repmat(sprintf("Bin%i", ii), length(binnedDatads{1, ii}), 1);
    if isempty(binvector)
        continue
    else
        if ii == 1
            DFFdscsvdata.Bin(1:DFFdsbineachlen(ii)) = binvector;
            DFFdscsvdata.DFF(1:DFFdsbineachlen(ii)) = binnedDatads{1, ii}.';
        elseif ii > 1
            DFFdscsvdata.Bin(DFFdsbineachlen(ii-1)+1:DFFdsbineachlen(ii)) = binvector;
            DFFdscsvdata.DFF(DFFdsbineachlen(ii-1)+1:DFFdsbineachlen(ii)) = binnedDatads{1, ii}.';
        end
    end
end

% AUC and DENS (ds and no ds):
AUCDENSbineachlen = cumsum(cellfun(@length, binnedData(2, :)));
AUCDENSbinlen = AUCDENSbineachlen(end);
AUCDENScsvdata = table(strings(AUCDENSbinlen, 1),nan(AUCDENSbinlen, 1),nan(AUCDENSbinlen, 1),nan(AUCDENSbinlen, 1),nan(AUCDENSbinlen, 1), ...
    'VariableNames', ["Bin" "AUC" "DENS" "AUCds" "DENSds"]);

for ii = 1:binsize
    binvector = repmat(sprintf("Bin%i", ii), length(binnedData{2, ii}), 1);
    if isempty(binvector)
        continue
    else
        if ii == 1
            AUCDENScsvdata.Bin(1:AUCDENSbineachlen(ii)) = binvector;
            AUCDENScsvdata.AUC(1:AUCDENSbineachlen(ii)) = binnedData{2, ii}.';
            AUCDENScsvdata.DENS(1:AUCDENSbineachlen(ii)) = binnedData{3, ii}.';
            AUCDENScsvdata.AUCds(1:AUCDENSbineachlen(ii)) = binnedDatads{2, ii}.';
            AUCDENScsvdata.DENSds(1:AUCDENSbineachlen(ii)) = binnedDatads{3, ii}.';
        elseif ii > 1
            AUCDENScsvdata.Bin(AUCDENSbineachlen(ii-1)+1:AUCDENSbineachlen(ii)) = binvector;
            AUCDENScsvdata.AUC(AUCDENSbineachlen(ii-1)+1:AUCDENSbineachlen(ii)) = binnedData{2, ii}.';
            AUCDENScsvdata.DENS(AUCDENSbineachlen(ii-1)+1:AUCDENSbineachlen(ii)) = binnedData{3, ii}.';
            AUCDENScsvdata.AUCds(AUCDENSbineachlen(ii-1)+1:AUCDENSbineachlen(ii)) = binnedDatads{2, ii}.';
            AUCDENScsvdata.DENSds(AUCDENSbineachlen(ii-1)+1:AUCDENSbineachlen(ii)) = binnedDatads{3, ii}.';
        end
    end
end


% MAX (max amplitude): 

MAXbineachlen = cumsum(cellfun(@length, binnedData(4, :)));
MAXbinlen = MAXbineachlen(end);
MAXcsvdata = table(strings(MAXbinlen, 1),nan(MAXbinlen, 1), 'VariableNames', ["Bin" "Max"]);

for ii = 1:binsize
    binvector = repmat(sprintf("Bin%i", ii), length(binnedData{4, ii}), 1);
    if isempty(binvector)
        continue
    else
        if ii == 1
            MAXcsvdata.Bin(1:MAXbineachlen(ii)) = binvector;
            MAXcsvdata.Max(1:MAXbineachlen(ii)) = binnedData{4, ii}.';
            MAXcsvdata.Maxds(1:MAXbineachlen(ii)) = binnedDatads{4, ii}.';
        elseif ii > 1
            MAXcsvdata.Bin(MAXbineachlen(ii-1)+1:MAXbineachlen(ii)) = binvector;
            MAXcsvdata.Max(MAXbineachlen(ii-1)+1:MAXbineachlen(ii)) = binnedData{4, ii}.';
            MAXcsvdata.Maxds(MAXbineachlen(ii-1)+1:MAXbineachlen(ii)) = binnedDatads{4, ii}.';
        end
    end
end

% Mean data: 

meantable = table(strings(binsize, 1), nan(binsize, 1), nan(binsize, 1), nan(binsize, 1), ...
    nan(binsize, 1), nan(binsize, 1), nan(binsize, 1), nan(binsize, 1), nan(binsize, 1), ...
    'VariableNames', ["Bin", "DFF", "DFFds", "AUC", "AUCds", "DENS", "DENSds", "Max", "Maxds"]);
binnames = strings(binsize, 1);
for ii = 1:binsize
    binnames(ii, 1) = sprintf("Bin%i", ii);
end
meantable.Bin = binnames;
meantable.DFF = meanDFF; meantable.DFFds = meanDFFds;
meantable.AUC = meanAUC; meantable.AUCds = meanAUCds;
meantable.DENS = meanDENS; meantable.DENSds = meanDENSds;
meantable.Max = meanMAX;meantable.Maxds = meanMAXds;

%% 6. Save all data in a struct variable: 

alldata.meantable = meantable;
alldata.DFFcsvdata = DFFcsvdata;
alldata.DFFdscsvdata = DFFdscsvdata;
alldata.MAXcsvdata = MAXcsvdata;
alldata.AUCDENScsvdata = AUCDENScsvdata;
alldata.binnedData = binnedData; % First row is DFF data, second row is AUC, third row is DENS, fourth row is max amplitude
alldata.binnedDatads = binnedDatads;
alldata.binbouts = binbouts;
alldata.params.threshold = threshold;
alldata.params.binsize  = binsize;
alldata.params.bins = bins;
alldata.params.binidx = binidx; 

%% 7. Save csv data in a specific folder: 

DFFname = strcat(savepath, "/binnedDFFdata.xlsx");
DFFdsname = strcat(savepath, "/binnedDFFdsdata.xlsx");
AUCDENSname = strcat(savepath, "/binnedAUCDENSdata.xlsx");
maxname = strcat(savepath, "/binnedmaxdata.xlsx");
meanname = strcat(savepath, "/binnedmeandata.xlsx");

writetable(DFFcsvdata,DFFname)
writetable(DFFdscsvdata,DFFdsname)
writetable(AUCDENScsvdata,AUCDENSname)
writetable(MAXcsvdata,maxname)
writetable(meantable,meanname)

end