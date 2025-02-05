function results = getripplesAlb(TDTdata,results,mypath)
% 
% % Prepare variables for data analysis:
% % id = '';
% % expcond = '';
% % prepost = '';
% thetaprop = [];
% thetaevpermin = [];
% LIAprop = [];
% LIAevpermin = [];
% meanrippleqt = [];
% totalrippleqt = [];
% meanLIAdur = [];
% totalLIAdur = [];
% meanripplefq_table = [];
% meanripdur_table = [];
% meansdpowrip_table = [];
% meanmaxpowrip_table = [];
% % mypath = '';
% 
% measurenames = {'Theta Time(%)' 'Theta Event per min(n/min)' 'LIA Time (%)' 'LIA Event per min (n/min)' 'mean rippleqt (n/LIA event)' 'total ripple qt (n)' 'mean LIA event duration (s)' 'total LIA event duration (s)' 'Mean Ripple Fq (Hz)' 'Mean Ripple Duration (ms)' 'Mean Power (SD)' 'Max Power (SD)' 'mypath'};
% Data = {thetaprop, thetaevpermin, LIAprop, LIAevpermin, meanrippleqt, totalrippleqt, meanLIAdur, totalLIAdur, meanripplefq_table ,meanripdur_table, meansdpowrip_table,meanmaxpowrip_table,mypath};
% 
% 
% 


    %% 1. Load saved file and other variables:

    % 1.1. Load variables:

    % Channels and ripple channel data:
    rippledata = results.lfp.datosfiltrados.';
    
    % LFP fs and PS fs:
    fs = TDTdata.streams.Wav1.fs;
    psfs = 2;
    
    % Frequencies:
    ripplefq = results.lfp.psprofilebueno.params.hfo_cf; %[90, 250]

    % Power spectrum of theta and ripples in both channels:
    
    rippleps = results.lfp.psprofilebueno.pS;

    % LIA threshold to detect ripples:
    LIAthresh = results.lfp.thetaanalisis.ratio.params.LIA;

    % Theta/delta ratios for each channel: 
    ratio = results.lfp.thetaanalisis.ratio.thetadelta;

    %Ripple duration y más cositas
    rippleduration = [0.015 0.25];
    noisethresh = 10;
    bandpassrg = [90 200];
    windur = 0.4;
    fqlimits = [30 250];
    

  

    % Generate time variable for LFP: 
    
    ts = (0:length(rippledata)-1)./fs;
    psts = (0:length(rippleps)-1)./psfs;
   

    %% 2. Get ripple events: 
    
   
    LIAidx = ratio < LIAthresh;
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

    boutidx_log = false(1, size(rippledata, 1));

    % Set the corresponding indices to true for each event
    for ii = 1:size(boutidx, 1)
        boutidx_log(boutidx(ii, 1):boutidx(ii, 2)) = true;
    end

    boutidx_log = boutidx_log.';  

    %% 3. Apply bandpass filter to raw signal - usar filtroFIR1.m (de Alberto): 

    pass1 = 90; % Hz
    pass2 = 200; % Hz
    
    [signRip_ripples, signEnv_ripples] = filtroFIR1(fs, double(rippledata), pass1, pass2);
    
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

           [cfs(:, :, ii),frq] = cwt(rippledata(ripplebouts(ii, 1)-lendif:ripplebouts(ii, 2)+lendif),fs, 'FrequencyLimits', fqlimits);
           tms = ((0:numel(rippledata(ripplebouts(ii, 1)-lendif:ripplebouts(ii, 2)+lendif))-1)/fs)-0.2;
       end

       idxtms = find(tms > 0,1,'first');
       for ii = 1:size(cfs, 3)
           [n, idx] = max(abs(cfs(:, idxtms, ii)));
           ripfrq(ii,1) = frq(idx);
           ripmaxpowtofrq(ii,1) = pow2db(n);
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
           disp('No hay ripples!!!!!!')
           results.rippleevent = NaN; %no hay. pelota
       else
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

       fig1name = strcat(mypath, '/ripple_events_cfs.jpg');
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

       % 6.3. Mean and max SD:

       for ii = 1:size(ripplebouts, 1)
           tempenvrip =    signEnv_ripples(ripplebouts(ii, 1):ripplebouts(ii, 2));
           ripmeansd(ii, 1) = mean((tempenvrip-envmean)/envsd); % mean sd of all sd values above threshold
           ripmaxsd(ii, 1) = (max(tempenvrip) - envmean)/envsd; % sd of max of envelope
       end

       % 6.4. Create data table with ripple measures:

       ripplenr = (1:size(ripplebouts, 1)).';

       % mouseid = strings(size(meanrippower, 1), 1);
       % expcond = zeros(size(meanrippower, 1), 1);
       event = zeros(size(meanrippower, 1), 1);
       durinsec = zeros(size(meanrippower, 1), 1);
       durinms = zeros(size(meanrippower, 1), 1);
       maxpow = zeros(size(meanrippower, 1), 1);
       meanpow = zeros(size(meanrippower, 1), 1);
       freq = zeros(size(meanrippower, 1), 1);
       record = zeros(size(meanrippower,1), 1);

       rippledatatable = table(event, durinsec, durinms, maxpow, meanpow, freq, ripmeansd, ripmaxsd,record, 'VariableNames', ["event", "durinsec", "durinms", "maxpow", "meanpow", "freq", "ripmeansd", "ripmaxsd", "record"]);
       % rippledata.expcond = num2str(rippledata.expcond);
       % rippledata.id = categorical(repmat(cellstr(results.id.mouseid), length(event), 1)); rippledata.expcond = categorical(repmat(cellstr(results.id.expcond), length(event), 1)); rippledata.event = ripplenr;
       rippledatatable.durinsec = rippledur.durinsec; rippledatatable.durinms = rippledur.durinms;
       rippledatatable.maxpow = ripmaxpowtofrq; rippledatatable.meanpow = meanrippower;
       % rippledata.record = repelem(nn, size(meanrippower,1)).';
       rippledatatable.freq = ripfrq; rippledatatable.ripmeansd = ripmeansd; rippledatatable.ripmaxsd = ripmaxsd;


       %% 7. Save data:

   % id = results.id.mouseid;
   % expcond = results.id.expcond;
   % prepost = results.id.time;
   % mypath = results.id.filename;
   totalrippleqt = size(event, 1);
   ripplerate = totalrippleqt/totalLIAbouts;
   meanripplefq_table = mean(ripfrq);
   meanripdur_table = mean(rippledur.durinms);
   meanpowrip_table = mean(meanrippower);
   meanmaxpowrip_table = mean(ripmaxpowtofrq);
   meansdperrec = mean(ripmeansd);
   maxsdperrec = mean(ripmaxsd);


   %% Theta/LIA proportion of time:
   
   Data = {totalrippleqt,ripplerate, meanripplefq_table ,meanripdur_table, meanpowrip_table,meanmaxpowrip_table,meansdperrec,maxsdperrec,mypath};
   measurenames = {'rippleqt' 'ripplerate' 'ripplefq' 'ripduration' 'meanpower' 'maxpower' 'meansd' 'maxsd' 'path'};
   results.rippleevent.ripplebouts = ripplebouts;
   results.rippleevent.rippledata = rippledatatable;
   %results.rippleevent.signEnv_complete = signEnv_complete;
   results.rippleevent.signEnv_ripples = signEnv_ripples;
   %results.rippleevent.signRip_complete = signRip_complete;
   results.rippleevent.signRip_ripples = signRip_ripples;
   results.LIAparams.ripplethresh = ripplethresh;
   results.rippleevent.cfs = cfs;
   results.rippleevent.cfsmeanrip = meanrip;
   results.rippleevent.cfsmeanripsm = meanripsm;
   results.rippleevent.params.fqlimits = fqlimits;
   results.rippleevent.frq = frq;

   eventsname = strcat(mypath, '/ripple_events.xlsx');
   writetable(rippledatatable, eventsname, 'WriteVariableNames', true)

   fig2name = strcat(mypath, '/mean_ripple_events.jpeg');
   saveas(fig2, fig2name);


   resultsname = char(strcat(mypath, '/Results.mat'));
   save(resultsname, 'results');
   
   close all
%   clearvars -except Data lfp_paths measurenames allripples bandpassrg noisethresh rippleduration windur nn fqlimits

   %% 8. Save csv files: 
   
       Data = cell2table(Data, 'VariableNames', measurenames);

       savename = mypath;
       dataname = strcat(savename, '/summaryrippledata.xlsx');
       rippledataname = strcat(savename, '/allrippleevents.xlsx');

       writetable(Data, dataname, 'WriteVariableNames', true);
       writetable(rippledatatable, rippledataname, 'WriteVariableNames', true);
       end
end % de la función