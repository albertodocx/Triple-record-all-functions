function Data = ICL_getripples(folderpath, varargin)

% INPUTS:
%  folderpath: main folder where data is stored in subfolders
%   OPTIONAL:
%  bandpassrg: matrix with the desired range where to filter for ripples. 
%   Default is 90 to 200 Hz
%  noisethresh: scalar number used to establish extreme hfo noise power 
%   values that are not to be included in the analysis. Default is 10
%   (Ej: 10*std de la envolvente)
%
% OUTPUTS: 
%  Data: ripple quantitative data 

p = inputParser;
addParameter(p,'bandpassrg', [90 200] ,@ismatrix);
addParameter(p, 'noisethresh', 10, @isnumeric);
addParameter(p, 'rippleduration', [0.015 0.25], @ismatrix);
addParameter(p, 'windur', 0.4, @isnumeric);
addParameter(p, 'fqlimits', [90 200], @ismatrix);
parse(p,varargin{:});
bandpassrg = p.Results.bandpassrg;
noisethresh = p.Results.noisethresh;
rippleduration = p.Results.rippleduration;
windur = p.Results.windur;
fqlimits = p.Results.fqlimits;


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
    if contains(listOfFolderNames(i), 'recording') && ~contains(listOfFolderNames(i), 'events') && ~contains(listOfFolderNames(i), 'continuous')
        lfp_paths = [lfp_paths listOfFolderNames(i)];
    else
        continue
    end
end

% Prepare variables for data analysis:
id = '';expcond = '';prepost = '';
totalLIAbouts = [];totalrippleqt = [];
ripplerate = [];meanripplefq_table = [];meanripdur_table = [];
meanpowrip_table = []; meansdperrec = [];
meanmaxpowrip_table = []; maxsdperrec = [];
path = '';

measurenames = {'Id' 'expcond' 'prepost'...
    'TotalLIAbouts' 'total ripple qt (n)'  'Mean Ripple Rate (Hz)' 'Mean Ripple Frequency (Hz)' 'Mean Ripple Duration (ms)' ...
    'Mean Power (dB)' 'Mean Power (SD)' 'Max Power (dB)' 'Max Power (SD)' 'path'};
Data = {id, expcond, prepost,...
       totalLIAbouts, totalrippleqt,ripplerate, meanripplefq_table ,meanripdur_table,...
       meanpowrip_table,meansdperrec,meanmaxpowrip_table,maxsdperrec,  path};



for nn = 1:length(lfp_paths)


    %% 1. Load saved file and other variables:

    filename = char(strcat(lfp_paths(nn), '\LFPdsdata.mat'));
    message = sprintf(['Processing file %i out of %i' ...
        ' ' ...
        ], nn, length(lfp_paths));
    disp(message)
    load(filename, 'results')

   savepath = results.id.filename;

    % 1.1. Load variables:

    % Channels and ripple channel data:
    ripplech = results.params.ripplech;
    thetach = results.params.thetach;
    ripchdata = results.LFPds.dsdata(:, ripplech);
    LFPdata = results.LFPds.dsdata;
    fs = results.LFPds.fs;

    
    % LFP fs and PS fs:
    lfpfs = results.LFPds.fs;
    psfs = 2;
    
    % Frequencies:
    thetafq = results.psProfile.params.theta_cf;
    ripplefq = results.psProfile.params.hfo_cf;

    % Power spectrum of theta and ripples in both channels:
    spthetaps = results.psProfile.pS(thetafq, :, ripplech);
    sprippleps = results.psProfile.pS(ripplefq, :, ripplech);

    slmthetaps = results.psProfile.pS(thetafq, :, thetach);
    slmrippleps = results.psProfile.pS(ripplefq, :, thetach);

    % LIA threshold to detect ripples:
    LIAthresh = results.thetaeventzscore.params.LIAthresh;

    % Theta/delta ratios for each channel: 
    spratio = results.thetaeventzscore.spthetafilt;
    slmratio = results.thetaeventzscore.slmthetafilt;

    % Generate time variable for LFP: 
    
    ts = (0:length(ripchdata)-1)./lfpfs;
    psts = (0:length(sprippleps)-1)./psfs;
    
    % Choose same ratio as the one used for theta event:
    
    eventch = results.thetaeventzscore.params.eventch;
    

    %% 2. Get ripple events: 
    
    if eventch == ripplech
        LIAidx = spratio < LIAthresh;
    elseif eventch == thetach
        LIAidx = slmratio < LIAthresh;
    end       
    
    [~, locs] = findpeaks(double(LIAidx));
    
    if LIAidx(1) == 1
        locs = [1, locs];
    end

           % Get onset and offset of each event saved in locs:

       for kk = 1:length(locs)
           if locs(kk) == locs(end)
               onset(kk) = locs(kk);
               if LIAidx(end) == 0
                   offset(kk) = locs(kk) + find(LIAidx(locs(kk):end) == 0, 1, 'first') - 1 ;
               else
                   offset(kk) = length(LIAidx);
               end
           else
               onset(kk) = locs(kk);
               offset(kk) = locs(kk) + find(LIAidx(locs(kk):end) == 0, 1, 'first') - 1 ;
           end
       end
    
    LIAbouts = [onset.', offset.'];

    LIAboutsts = psts(LIAbouts);
    
    totalLIAbouts = sum(LIAboutsts(:, 2)-LIAboutsts(:,1)); % seg

    for ii = 1:size(LIAboutsts, 1)
        boutidx(ii, 1) = find(ts >= LIAboutsts(ii, 1),1, 'first');
        boutidx(ii, 2) = find(ts <= LIAboutsts(ii, 2),1, 'last');
    end

    % Convert bouts into logical time series: 

    boutidx_log = false(1, size(LFPdata, 1));

    % Set the corresponding indices to true for each event
    for ii = 1:size(boutidx, 1)
        boutidx_log(boutidx(ii, 1):boutidx(ii, 2)) = true;
    end

    boutidx_log = boutidx_log.';    


    %% 3. Apply bandpass filter to raw signal - usar filtroFIR1.m (de Alberto): 

    % Of the whole channel and LFPdata: 

    pass1 = bandpassrg(1); 
    pass2 = bandpassrg(2); 
    
    [signRip_complete, signEnv_complete] = filtroFIR1(lfpfs, LFPdata, pass1, pass2);
    [signRip_ripples, signEnv_ripples] = filtroFIR1(lfpfs, ripchdata, pass1, pass2);
   
    % Control de valores extremos de envolvente: 

    hfonoise = mean(signEnv_ripples + noisethresh*std(signEnv_ripples));
    signEnv_ripples(signEnv_ripples > hfonoise) = nan;
    envmean = mean(signEnv_ripples, 'omitnan');
    envsd = std(signEnv_ripples, 'omitnan');
    ripplethresh = mean(signEnv_ripples, 'omitnan') + 3*std(signEnv_ripples, 'omitnan');


    %% 4. Isolate ripple bouts:

    % Find ripple bouts: 
    rippleidx = signEnv_ripples > ripplethresh & boutidx_log;
    [~, locs] = findpeaks(double(rippleidx));
    ripplebouts = zeros(size(locs, 1), 2);

      for kk = 1:length(locs)
           if locs(kk) == locs(end)
               ripplebouts(kk, 1) = locs(kk);
               if rippleidx(end) == 0
                   ripplebouts(kk, 2) = locs(kk) + find(rippleidx(locs(kk):end) == 0, 1, 'first') - 1 ;
               else
                   ripplebouts(kk, 2) = length(rippleidx);
               end
           else
               ripplebouts(kk, 1) = locs(kk);
               ripplebouts(kk, 2) = locs(kk) + find(rippleidx(locs(kk):end) == 0, 1, 'first') - 1 ;
           end
       end
    
       if isempty(ripplebouts)
           clearvars -except Data lfp_paths measurenames allripples bandpassrg noisethresh rippleduration windur nn
           continue
       end
       rippledur.durinsec = (ripplebouts(:, 2) - ripplebouts(:, 1))./fs;
       rippledur.durinms = rippledur.durinsec*1000;


       %% 5. Discarding criteria:

       % Minimum ripple duration 15-250 ms (Tingley, D., McClain, K., Kaya, E.,
       % Carpenter, J., & Buzsáki, G. (2021).
       % A metabolic function of the hippocampal sharp wave-ripple. Nature, 597(7874), 82-86.)

       mindur = rippleduration(1); maxdur = rippleduration(2);
       mindurindx = rippledur.durinsec >= mindur & rippledur.durinsec <= maxdur;
       ripplebouts = ripplebouts(mindurindx, :);
       rippledur.durinsec = rippledur.durinsec(mindurindx);
       rippledur.durinms = rippledur.durinms(mindurindx);

       %% 6. Calculate wavelet-Morlet spectrum of each ripple event:

       winlength = windur*fs; % 400 ms de ventana de duración de los ripples

       % tiledlayout(figsize,figsize)
       for ii = 1:size(ripplebouts, 1)
           lencfs = ripplebouts(ii, 2)-ripplebouts(ii, 1);
           lendif = (winlength - lencfs)/2;
           if ripplebouts(ii, 1) - lendif < 0
               continue
           end
           [cfs(:, :, ii),frq] = cwt(LFPdata(ripplebouts(ii, 1)-lendif:ripplebouts(ii, 2)+lendif, 3),fs, 'FrequencyLimits', fqlimits);
           tms = ((0:numel(LFPdata(ripplebouts(ii, 1)-lendif:ripplebouts(ii, 2)+lendif, 3))-1)/fs)-0.2;

           %     nexttile
           %     surface(tms, frq,abs(cfs(:, :, ii)))
           %     title(ii)
           %     axis tight
           %     shading flat
       end

       % 6.1. Get intrinsic frequency and max power when time = 0 (center of the
       % ripple)

       for ii = 1:size(cfs, 3)
           [n, idx] = max(abs(cfs(:, tms == 0, ii)));
           ripfrq(ii, 1) = frq(idx);
           ripmaxpowtofrq(ii, 1) = pow2db(n);
       end

       fqthresh = (ripfrq > pass1 & ripfrq < pass2);
       maxpowthresh = (ripmaxpowtofrq >= (mean(ripmaxpowtofrq)-3*std(ripmaxpowtofrq)) & ripmaxpowtofrq <= (3*std(ripmaxpowtofrq)+mean(ripmaxpowtofrq)));

       discardripevents = (fqthresh + maxpowthresh) == 2;
       ripplebouts = ripplebouts(discardripevents, :);
       cfs = cfs(:, :, discardripevents);
       ripfrq = ripfrq(discardripevents);
       ripmaxpowtofrq = ripmaxpowtofrq(discardripevents);
       rippledur.durinsec = rippledur.durinsec(discardripevents);
       rippledur.durinms = rippledur.durinms(discardripevents);

       if isempty(ripplebouts)
           clearvars -except Data lfp_paths measurenames allripples bandpassrg noisethresh rippleduration windur nn
           continue
       end

       % Check ripples:

       figsize = floor(sqrt(size(ripplebouts, 1)));
       if figsize > 10
           figsize = 10;
       end
       tiledlayout(figsize, figsize);
       for ii = 1:(figsize*figsize)
           nexttile
           surface(tms, frq,abs(cfs(:, :, ii)))
           title(ii)
           axis tight
           shading flat
       end

       fig1 = figure(1);
       newWidth = 1600; 
       newHeight = 900; 
       set(fig1, 'Position', [100, 100, newWidth, newHeight])

       fig1name = strcat(savepath, '\ripple_events_cfs.jpeg');
       saveas(fig1, fig1name);


       % Get mean power of this recording:

       meanrip = zeros(size(frq, 1), size(cfs, 2));
       for ii = 1:length(frq)
           for jj = 1:size(cfs, 3)
               temp(jj, :) = cfs(ii, :, jj);
           end
           meanrip(ii, :) = mean(pow2db(abs(temp)), 1);
       end

       meanripsm = zeros(size(meanrip, 1), size(meanrip, 2));
       for ii = 1:size(meanrip, 2)
           meanripsm(:, ii) = smooth(meanrip(:, ii), 'moving',500);
       end

       tscfs = linspace(-0.2, 0.2, size(cfs, 2));

       fig2 = figure(2);
       subplot(2,1,1)
       surface(tscfs,frq,meanrip)
       axis tight
       shading flat
       xlabel("Time (s)")
       ylabel("Frequency (Hz)")
       set(gca,"yscale","log")

       subplot(2,1,2)
       surface(tscfs,frq,meanripsm)
       %contourf(tms, frq, pow2db(abs(cfs)))
       axis tight
       shading flat
       xlabel("Time (s)")
       ylabel("Frequency (Hz)")
       set(gca,"yscale","log")

       % 6.2. Mean ripple power:

       for ii = 1:size(cfs, 3)
           onset = find(tms >= 0-(rippledur.durinsec(ii)/2),1, 'first');
           offset = find(tms <= 0+(rippledur.durinsec(ii)/2),1, 'last');
           x = onset:offset;
           y = find(frq == ripfrq(ii));
           if y < length(frq) && y > 1
               meanrippower(ii, 1) = mean(mean(pow2db(abs(cfs(y-1:y+1, x, ii))), 2), 1);
           elseif y == length(frq)
               meanrippower(ii, 1) = mean(mean(pow2db(abs(cfs(y-2:y, x, ii))), 2), 1);
           elseif y == 1
               meanrippower(ii, 1) = mean(mean(pow2db(abs(cfs(y:y+2, x, ii))), 2), 1);
           end
       end


       % 6.3. Entropy:

       % Normalized power spectrum:
       %
       % for ii = 1:size(cfs, 3)
       %
       % end

       % 6.4. Mean and max SD:

       for ii = 1:size(ripplebouts, 1)
           tempenvrip =    signEnv_ripples(ripplebouts(ii, 1):ripplebouts(ii, 2));
           ripmeansd(ii, 1) = mean((tempenvrip-envmean)/envsd); % mean sd of all sd values above threshold
           ripmaxsd(ii, 1) = (max(tempenvrip) - envmean)/envsd; % sd of max of envelope
       end

       % 6.5. Create data table with ripple measures:

       ripplenr = [1:size(ripplebouts, 1)].';

       mouseid = strings(size(meanrippower, 1), 1);
       expcond = zeros(size(meanrippower, 1), 1);
       event = zeros(size(meanrippower, 1), 1);
       durinsec = zeros(size(meanrippower, 1), 1);
       durinms = zeros(size(meanrippower, 1), 1);
       maxpow = zeros(size(meanrippower, 1), 1);
       meanpow = zeros(size(meanrippower, 1), 1);
       freq = zeros(size(meanrippower, 1), 1);
       record = zeros(size(meanrippower,1), 1);

       rippledata = table(mouseid, expcond, event, durinsec, durinms, maxpow, meanpow, freq, ripmeansd, ripmaxsd,record, 'VariableNames', ["id", "expcond", "event", "durinsec", "durinms", "maxpow", "meanpow", "freq", "ripmeansd", "ripmaxsd", "record"]);
       rippledata.expcond = num2str(rippledata.expcond);
       rippledata.id = categorical(repmat(cellstr(results.id.mouseid), length(event), 1)); rippledata.expcond = categorical(repmat(cellstr(results.id.expcond), length(event), 1)); rippledata.event = ripplenr;
       rippledata.durinsec = rippledur.durinsec; rippledata.durinms = rippledur.durinms;
       rippledata.maxpow = ripmaxpowtofrq; rippledata.meanpow = meanrippower;
       rippledata.record = repelem(nn, size(meanrippower,1)).';
       rippledata.freq = ripfrq; rippledata.ripmeansd = ripmeansd; rippledata.ripmaxsd = ripmaxsd;


       %% 7. Save data:

   id = results.id.mouseid;
   expcond = results.id.expcond;
   prepost = results.id.time;
   path = results.id.filename;
   totalrippleqt = size(event, 1);
   ripplerate = totalrippleqt/totalLIAbouts;
   meanripplefq_table = mean(ripfrq);
   meanripdur_table = mean(rippledur.durinms);
   meanpowrip_table = mean(meanrippower);
   meanmaxpowrip_table = mean(ripmaxpowtofrq);
   meansdperrec = mean(ripmeansd);
   maxsdperrec = mean(ripmaxsd);


   %% Theta/LIA proportion of time:
   
   Data(nn, :) = {id, expcond, prepost,totalLIAbouts, totalrippleqt,ripplerate, meanripplefq_table ,meanripdur_table, meanpowrip_table,meanmaxpowrip_table,meansdperrec,maxsdperrec,  path};

   results.rippleevent.ripplebouts = ripplebouts;
   results.rippleevent.rippledata = rippledata;
   results.rippleevent.signEnv_complete = signEnv_complete;
   results.rippleevent.signEnv_ripples = signEnv_ripples;
   results.rippleevent.signRip_complete = signRip_complete;
   results.rippleevent.signRip_ripples = signRip_ripples;
   results.LIAparams.ripplethresh = ripplethresh;
   results.rippleevent.cfs = cfs;
   results.rippleevent.cfsmeanrip = meanrip;
   results.rippleevent.cfsmeanripsm = meanripsm;
   results.rippleevent.params.fqlimits = fqlimits;

   eventsname = strcat(savepath, '\ripple_events.xlsx');
   writetable(rippledata, eventsname, 'WriteVariableNames', true)

   fig2name = strcat(savepath, '\mean_ripple_events.jpeg');
   saveas(fig2, fig2name);

    if nn == 1
        allripples = rippledata;
    elseif nn > 1 && ~exist('allripples', 'var') == 1
        allripples = rippledata;
    else 
        allripples = [allripples; rippledata];
    end

    allripples.id = categorical(cellstr(allripples.id));

   resultsname = char(strcat(lfp_paths(nn), '\LFPdsdata.mat'));
   save(resultsname, 'results');
   
   close all
   clearvars -except Data lfp_paths measurenames allripples bandpassrg noisethresh rippleduration windur nn
   


end % loop each recording

   %% 8. Save csv files: 
   
       Data = cell2table(Data, 'VariableNames', measurenames);

       savename = uigetdir();
       dataname = strcat(savename, '\summaryrippledata.xlsx');
       allripplename = strcat(savename, '\allrippleevents.xlsx');

       writetable(Data, dataname, 'WriteVariableNames', true);
       writetable(allripples, allripplename, 'WriteVariableNames', true);



end % function
