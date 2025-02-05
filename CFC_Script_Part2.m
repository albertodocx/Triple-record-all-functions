
%%%%%%%%%%%%%%%%    FIBER PHOT. vs CONTINUOUS VARIABLES     %%%%%%%%%%%%%%%

%% 1. Extract all files in batch: 

% Specify the root folder
startFolder = uigetdir();

% Use the dir function to get information about all files in the root folder and its subfolders
fileList = dir(fullfile(startFolder, '**', 'results.mat'));

% Display the list of matching files
mypaths = strings(size(fileList, 1), 1);
for ii = 1:length(fileList)
    % Check if the entry is a file (not a directory)
    if ~fileList(ii).isdir
        mypaths(ii) = convertCharsToStrings((fullfile(fileList(ii).folder, fileList(ii).name)));
    end
end

%% %%%%%%%%%%%%%%%    BINNING DATA    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2.1. Check and determine continuous data thresholds: 

% Behavioral variables should be stored in the same structure
% (results.Behavior): 
behfs  = 15;
% 5min:
requiredlength = round(5*60*behfs); % 5min * fs
vel = nan(length(mypaths), requiredlength);
accel = nan(length(mypaths), requiredlength);
MDB = nan(length(mypaths), requiredlength);

for ii = 1:length(mypaths)
    load(mypaths(ii))
    vel(ii, 1:length(results.Behavior.Velocity)) = results.Behavior.Velocity;
    accel(ii, 1:length(results.Behavior.Velocity)) = [0; CFC_calcaccel(results.Behavior.Velocity)];
    MDB(ii, 1:size(results.Behavior.Location, 1)) = CFC_calcMDB(results.Behavior.Location, behfs);
end

% Find bin-thresholds, example for velocity: 
velsm = arrayfun(@(row) smooth(vel(row, :), behfs), 1:size(vel, 1), 'UniformOutput', false);
velsm = cell2mat(velsm).';velmax = max(velsm,[], 1); velmin = min(velsm, [], 1);
velmax(velmax < 0) = nan; velmin(velmin<0) = nan;
[a, ~] = mle(velmax(~(isnan(velmax))),'Distribution' ,'InverseGaussian');
[b, ~] = mle(velmin(~(isnan(velmin))),'Distribution' ,'InverseGaussian'); 
mumax = a(1);lambdamax = a(2); alphamax = 0.05;
mumin = b(1);lambdamin = b(2); alphamin = 0.05;
velthresh(2) = icdf('InverseGaussian', 1 - alphamax, mumax, lambdamax);
velthresh(1) = icdf('InverseGaussian', 1 - alphamin, mumin, lambdamin);

% Acceleration: 
accelsm = arrayfun(@(row) smooth(accel(row, :), behfs), 1:size(accel, 1), 'UniformOutput', false);
accelsm = cell2mat(accelsm).'; accelmean = mean(accelsm, 1, 'omitnan');
accelmin = min(accelsm, [], 1);accelmax = max(accelsm, [], 1);
accelthresh(1) = mean(accelmin) - 1.96*std(accelmin);
accelthresh(2) = mean(accelmax) + 1.96*std(accelmax);

% plot(accelmin)
% yline(mean(accelmin)-1.96*std(accelmin))

% figure(1)
% set(gcf,'Position',[100 50 1000 1000])
% tiledlayout(6, 4)
% for ii = 1:22
%     nexttile
%     plot(velsm(ii, :))
%     yline(value)
%     ylim([0 20])
% end
% 

% 3 min: 
% requiredlength3 = round(3*60*behfs); % 3min * fs. 
% vel3 = vel(:, 1:requiredlength3);
% accel3 = accel(:, 1:requiredlength3);
% MDB3 = MDB(:, 1:requiredlength3);

% Calculate thresholds for the bins continuous variables: 

% % 5 min
% [vel_thresh(1), vel_thresh(2)] = BEH_findbinthresholds(vel, behfs);
% [acc_thresh(1), acc_thresh(2)] = BEH_findbinthresholds(accel, behfs);
% % Para MDB no tiene sentido ya que ya tenemos los límites: 0 a 12,5 cm 
% 
% % 3 min
% [vel_thresh3(1), vel_thresh3(2)] = BEH_findbinthresholds(vel3, behfs);
% [acc_thresh3(1), acc_thresh3(2)] = BEH_findbinthresholds(accel3, behfs);
% 


%% 2.2. Loop through all files to get bins for all results files: 

savepathfolders = strings(length(mypaths), 1);
for ii = 1:length(mypaths)
    savepathfolders(ii) = convertCharsToStrings(fileList(ii).folder);
end

% Establish constant variables for all files: 
binsize = 8;
binnedData = cell(1, length(mypaths));
binnedData_accel = cell(1, length(mypaths));
binnedData3 = cell(1, length(mypaths));
% Loop through all files: 
for ii = 1:length(mypaths)
    load(mypaths(ii))
    %%%%%%%%%%%%%%%% 5 MIN %%%%%%%%%%%%%%%%%%
    % FP Signals:
    DFF = results.FP.Signals.DFFModZscore;
    fs = results.FP.params.fs;
    % Behavior or continuous signal to bin (ejemplos):
    vel = results.Behavior.Velocity;
    behfs = results.Behavior.Fs;
    % Extract binned data:
    binnedData{ii} = BEH_PHOT_binvar(vel, DFF, behfs, fs, velthresh, binsize, savepathfolders(ii));
    %     binnedData_accel{ii} = BEH_PHOT_binvar(accel, DFF, behfs, fs, accelthresh, binsize, savepathfolders(ii));
    %%%%%%%%%%%%%%%% 3 MIN %%%%%%%%%%%%%%%%%%
    DFF3 = DFF(1:round(fs*60*3));
    vel3 = vel(1:round(behfs*60*3));
    binnedData3{ii} = BEH_PHOT_binvar(vel3, DFF3, behfs, fs, velthresh, binsize, savepathfolders(ii));   
end

% Join all tables with same info: 

allmean = binnedData{1, 1}.meantable;
allDFF = binnedData{1,1}.DFFcsvdata;
allDFFds = binnedData{1,1}.DFFdscsvdata;
allAUCDENS = binnedData{1,1}.AUCDENScsvdata;
allMAX = binnedData{1,1}.MAXcsvdata;

for ii = 2:length(mypaths)
    allmean = vertcat(allmean, binnedData{1, ii}.meantable);
    allDFF = vertcat(allDFF, binnedData{1, ii}.DFFcsvdata);    
    allDFFds = vertcat(allDFFds, binnedData{1, ii}.DFFdscsvdata);    
    allAUCDENS = vertcat(allAUCDENS, binnedData{1, ii}.AUCDENScsvdata);    
    allMAX = vertcat(allMAX, binnedData{1, ii}.MAXcsvdata);    
end


%% 2.3. Save data for stats :) 

alldatapath = uigetdir();
allmeanname = strcat(alldatapath, "/BINSallmeandata.xlsx");
allDFFname = strcat(alldatapath, "/BINSallDFFdata.xlsx");
allDFFdsname = strcat(alldatapath, "/BINSallDFFdsdata.xlsx");
allAUCDENSname = strcat(alldatapath, "/BINSallAUCDENSdata.xlsx");
allMAXname = strcat(alldatapath, "/BINSallMAXdata.xslx"); 

writetable(allmean, allmeanname)
writetable(allDFF, allDFFname)
writetable(allDFFds, allDFFdsname)
writetable(allAUCDENS, allAUCDENSname)
writetable(allMAX, allMAXname)


%% 2.4. Create PETHs of any variable: 
 % check for options to modify the PETH in the "optional inputs" section of the
 % function
% El max length de los registros por defecto se corresponde con el del CFC.
% Se puede modificar para otros registros más largos. Es simplemente un
% control de que no se excedan las muestras (cálculos y redondeos de fs)

fsvel = results.Behavior.Fs; % fs de var de interés/q en caso de CFC es igual que la velocidad y otras vars de conducta
vel = smooth(results.Behavior.Velocity, fscontvar, 'moving'); % extraer variable de interés
accel = smooth(diff(vel), fscontvar, 'moving');
% accel = diff(smooth(results.Behavior.Velocity), fscontvar, 'moving');
freezsamp = results.Behavior.Event.Time; % freezing 
fsev = fsvel;
savepath = uigetdir(); % ruta donde guardar los datos
savepathPETH = uigetdir(); % ruta donde guardar los datos
DFF = results.FP.Signals.DFFModZscore;
fsDFF = results.FP.params.fs;
shocks_evsamp = [round(fsev*3*60); round(fsev*3.5*60); round(fsev*4*60)];
shocks_evsamp(:, 2) = round(shocks_evsamp+fsev.*2); 

%%%%%%%%%%%%%%%%%%%% Examples: 
%  >> freezing events and DFF
results.Analysis.freezandDFF = BEH_PETHonset(DFF, fsDFF, freezsamp, fsev,  savepath);
%  >> movement events and DFF (para FP + LFP, cuando tengamos eventos de theta)
%  >> shocks and DFF
results.Analysis.shocksandDFF = BEH_PETHonset(DFF, fsDFF, shocks_evsamp, fsev,  savepath);
%  >> freezing events and velocity/acceleration
results.Analysis.freezandvel = BEH_PETHonset(vel, fsvel, freezsamp, fsev,  savepathPETH);
results.Analysis.freezandaccel = BEH_PETHonset(accel, fsvel, freezsamp, fsev,  savepathPETH);
%  >> shocks and velocity/acceleration
results.Analysis.shocksandvel = BEH_PETHonset(vel, fsvel, shocks_evsamp, fsev,  savepathPETH);
results.Analysis.shocksandvel = BEH_PETHonset(accel, fsvel, shocks_evsamp, fsev,  savepathPETH);


%% 2.5. SAME AS 2.4 BUT IN A FOR LOOP FOR ALL FILES AND MANY COMBINATIONS OF VARIABLES: 

% Shock samples: generates a matrix like events but with the timing of
% shocks and a duration of 2 seconds (30+ samples according to behfs)

% for ii = 1:length(mypaths)
velfreezname = 'velandfreez'; velshname = 'velandshocks';
accelfreezname = 'accelandfreez'; acshname = 'accelandshocks';
dfffreezname = 'dffandfreez'; dffshname = 'dffandshocks';
% mdbfreezname = 'mdbandfreez'; mdbshname = 'mdbandshocks';
for ii = 1:1
    load(mypaths(ii))
    %%%%%%%%%%%%%%%% 5 MIN %%%%%%%%%%%%%%%%%%
    % FP Signals:
    DFF = results.FP.Signals.DFFModZscore;
    fsDFF = results.FP.params.fs;
    % Behavior or continuous signal:
    fsvel = results.Behavior.Fs;
    vel = smooth(results.Behavior.Velocity, fsvel, 'moving');
    accel = diff(vel);
%     MDB = CFC_calcMDB(results.Behavior.Location, fsvel).';
    % EVENTS: 
    freezsamp = results.Behavior.Event.Time; % freezing
    fsev = fsvel;
    shocks_evsamp = [round(fsev*3*60); round(fsev*3.5*60); round(fsev*4*60)];
    shocks_evsamp(:, 2) = round(shocks_evsamp+fsev.*2);
    % shocks_evsamp: already defined before the loop as is constant, though
    % only in the 3rd day. 
    % SAVEPATH:
     savepath = results.FP.path;
    % EXTRACT PETHs: 
    results.Analysis.freezandDFF = BEH_PETHonset(DFF, fsDFF, freezsamp, fsev, dfffreezname, savepath);
    results.Analysis.freezandvel = BEH_PETHonset(vel, fsvel, freezsamp, fsev, velfreezname, savepath);
    results.Analysis.freezandaccel = BEH_PETHonset(accel, fsvel, freezsamp, fsev,accelfreezname,  savepath);
%     results.Analysis,freezandMDB = BEH_PETHonset(MDB, fsvel, freezsamp, fsev, mdbfreezname, savepath);
    results.Analysis.shocksandDFF = BEH_PETHonset(DFF, fsDFF, shocks_evsamp, fsev,dffshname,  savepath);
    results.Analysis.shocksandvel = BEH_PETHonset(vel, fsvel, shocks_evsamp, fsev,velshname,  savepath);
    results.Analysis.shocksandaccel = BEH_PETHonset(accel, fsvel, shocks_evsamp, fsev,acshname,  savepath);
%     results.Analysis.shocksandMDB = BEH_PETHonset(MDB, fsvel, shocks_evsamp, fsev, mdbshname, savepath);
    save(savepath, 'results')
end



%% 2.6. SAME AS 2.5 BUT NO NORMALIZATION: 

savepath = uigetdir();
contvar = results.Behavior.Velocity;
contvar = smooth(diff(contvar), fscontvar, 'moving');
combiname = 'accelandfreez';
mydata = BEH_PETHonset_realvalues(contvar, fscontvar, evsamp, fsev, combiname,  savepath);





