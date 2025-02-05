
function Data = ICL_thetaevent(folderpath, varargin)


% INPUTS:
%  folderpath: main folder where data is stored in subfolders
%   OPTIONAL:
%  thetafqrng: matrix with the desired theta range. Default is 5.5 to 8.5
%  Hz (igual que Alberto)
%  deltafqrng: matrix with the desired delta range. Default is 1 to 3.5 Hz
%  (igual que Alberto)
%  LIAthresh: ripples are detected within periods under this threshold 
%  (normalized from 0 to 1).
%  runthresh: real run theta events are detected within periods above this
%  threshold (normalized from 0 to 1). 
%
% OUTPUTS: 
%  Data: a table containing quantitative data about theta events like mean
%  power of each frequency band (gamma and theta) and ID info of each file.


p = inputParser;
addParameter(p,'thetafqrng', [5.5 8.5] ,@ismatrix);
addParameter(p,'deltafqrng', [1 3.5] ,@ismatrix);
addParameter(p,'LIAthresh', 0.2 ,@isnumeric);
addParameter(p,'runthresh', 0.3 ,@isnumeric);
parse(p,varargin{:});
thetafqrng = p.Results.thetafqrng;
deltafqrng = p.Results.deltafqrng;
LIAthresh = p.Results.LIAthresh;
runthresh = p.Results.runthresh;



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

% 1.1. Prepare csv tables:

id = '';
expcond = '';
prepost = '';
eventnr = [];
evdurtot = [];
evdurmean = [];
totmean_theta = [];
totmean_gamma = [];
totmean_gammahigh = [];
totmean_gammalow = [];
path = '';

   measurenames = {'Id' 'expcond' 'prepost' 'eventnr' 'tot eventdur (s)' 'mean eventdur (s)' 'thetamean (dB)' 'gammamean (dB)' 'gammahighmean (dB)' 'gammalowmean(dB)' 'path'};
   Data = {id, expcond, prepost, eventnr, evdurtot, evdurmean ,totmean_theta, totmean_gamma,totmean_gammahigh,totmean_gammalow,  path};


for nn = 1:length(lfp_paths)


    %% 1. Load saved file:

    filename = char(strcat(lfp_paths(nn), '\LFPdsdata.mat'));
    message = sprintf(['Processing file %i out of %i' ...
        ' ' ...
        ], nn, length(lfp_paths));
    disp(message)
    load(filename, 'results')


    %% 2. Extract and generate main variables:

    % 2.1. SLM channel:

    thetach = results.params.thetach;
    thetach_profile = results.psProfile.pS(:, :, thetach);
    
    thetafq = (results.psProfile.pSfreqs >= thetafqrng(1) & results.psProfile.pSfreqs <= thetafqrng(2)); % fq range
    thetach_psth = thetach_profile(thetafq, :); 
    deltafq = (results.psProfile.pSfreqs >= deltafqrng(1) & results.psProfile.pSfreqs <= deltafqrng(2)); % based on Malezieux, M., Kees, A. L., & Mulle, C. (2020)
    thetach_psde= thetach_profile(deltafq, :);

    if size(thetach_psde, 1) > size(thetach_psth, 1)
        thetach_psde = thetach_psde(1:size(thetach_psth, 1), :);
    elseif size(thetach_psde, 1) < size(thetach_psth, 1)
        thetach_psth = thetach_psth(1:size(thetach_psde, 1), :);
    end
    
%     psfs = results.psProfile.params.tWin/1000;
% %     time_spec = [0:size(thetach_psth, 2)-1]./psfs;
%     time_spec = (0:size(thetach_psth, 2)-1)./2; % 2 muestras por seg

% 2.2. SP channel: 

    ripplech = results.params.ripplech;
    ripplech_profile = results.psProfile.pS(:, :, ripplech);
    
    sp_thetafq = (results.psProfile.pSfreqs >= thetafqrng(1) & results.psProfile.pSfreqs <= thetafqrng(2)); % fq range
    ripplech_psth = ripplech_profile(sp_thetafq, :); 
    sp_deltafq = (results.psProfile.pSfreqs >= deltafqrng(1) & results.psProfile.pSfreqs <= deltafqrng(2)); % based on Malezieux, M., Kees, A. L., & Mulle, C. (2020)
    ripplech_psde= ripplech_profile(sp_deltafq, :);

    if size(ripplech_psde, 1) > size(ripplech_psth, 1)
        ripplech_psde = ripplech_psde(1:size(ripplech_psth, 1), :);
    elseif size(ripplech_psde, 1) < size(ripplech_psth, 1)
        ripplech_psth = ripplech_psth(1:size(ripplech_psde, 1), :);
    end
    
    psfs = 2; 
    time_spec = (0:size(ripplech_psth, 2)-1)./2; % 2 muestras por seg
    freqs = results.psProfile.pSfreqs;

    
    %% 3. Ratio theta/delta:

    % 3.1. SLM channel:
    thetamean = mean(thetach_psth, 1); deltamean = mean(thetach_psde, 1);
    ratiotd = thetamean./deltamean;
    if sum(isinf(ratiotd)) > 0
        ratiotd(ratiotd == Inf) = nan;
    end

    % Código de Alberto para sgolayfilt:
    ratiofilt = movmean( movmean( sgolayfilt(ratiotd,50,51), 20), 20); 
    ratiofilt = rescale(ratiofilt, 0, 1);
    ratiosc = rescale(ratiotd, 0, 1); 

    % Apply sgolayfilt + rescale also to theta power (como Alberto):

    thetameansc = rescale(thetamean, 0, 1);
    thetafilt =  movmean( movmean( sgolayfilt(thetamean,50,51), 20), 20);
    thetafilt = rescale(thetafilt, 0, 1);


    % 3.2. SP channel:

    sp_thetamean = mean(ripplech_psth, 1); sp_deltamean = mean(ripplech_psde, 1);
    sp_ratiotd = sp_thetamean./sp_deltamean;
    if sum(isinf(sp_ratiotd)) > 0
        sp_ratiotd(sp_ratiotd == Inf) = nan;
    end

    sp_ratiosc = rescale(sp_ratiotd, 0, 1); % normalización como Alberto a 0-1.

    % Código de Alberto para sgolayfilt:
    sp_ratiofilt = movmean( movmean( sgolayfilt(sp_ratiotd,50,51), 20), 20);
    sp_ratiofilt = rescale(sp_ratiofilt, 0, 1);
    %     sp_ratiofiltmean = mean(sp_ratiofilt, 'omitnan');
    %     sp_ratiofiltstd = std(sp_ratiofilt);

    sp_thetameansc = rescale(sp_thetamean, 0, 1);
    sp_thetafilt =  movmean( movmean( sgolayfilt(sp_thetamean,50,51), 20), 20);
    sp_thetafilt = rescale(sp_thetafilt, 0, 1);


    %% 4. Apply criteria to detect run-theta events

    %%%%%%%%%%%%%%%% CRITERIA TO DETECT RUN-THETA EVENTS %%%%%%%%%%%%%%%%%%%%%
    % 1. Thresholds: 0.5 for run theta. This threshold is applied to filtered
    % ratio (env, sgolayfilt).
    % 2. Events that are less than 1 sec apart from each other are considered
    % part of the same event.
    % 3. Events shorter than 1 sec are not considered as theta events

    % 4.1.  THRESHOLDS COMO ALBERTO:

    % Compare SP and SLM to check which channel to choose for theta calc.:

    imrg = results.psProfile.pSfreqs >= 1 & results.psProfile.pSfreqs <= 15;
    imgfreqs = freqs(imrg);

    fig1 = figure(nn);

    tiledlayout(3, 2)
    % Power specturm SP and SLM:
    %SP
    nexttile
    imagesc(time_spec, imgfreqs, pow2db(ripplech_profile(imrg, :)))
    title('SP - Power Spectrum (dB)')
    ylabel('Fq (Hz)')
    xlabel('Time (s)')
    colorbar
    colormap('jet')
    set(gca, 'YDir', 'normal')
    clim([0 100])
    % SLM: 
    nexttile
    imagesc(time_spec, imgfreqs, pow2db(thetach_profile(imrg, :)))
    title('SLM - Power Spectrum (dB)')
    ylabel('Fq (Hz)')
    xlabel('Time (s)')
    colorbar
    colormap('jet')
    set(gca, 'YDir', 'normal')
    clim([0 100])


    % SP
    nexttile % Plot normalized theta power
    plot(time_spec, sp_thetameansc, 'Color' ,'#dad0d8')
    hold on
    plot(time_spec, sp_thetafilt, 'blue')
    title('SP. Normalized Theta Power')
    xlabel('Time (s)')
    hold off
    % SLM
    nexttile % Plot normalized theta power
    plot(time_spec, thetameansc, 'Color' ,'#dad0d8')
    hold on
    plot(time_spec, thetafilt, 'blue')
    title('SLM. Normalized Theta Power')
    xlabel('Time (s)')
    hold off

    % SP:
    nexttile % Plot normalized ratio power
    plot(time_spec, sp_ratiosc, 'Color' ,'#dad0d8')
    hold on
    plot(time_spec, sp_ratiofilt, 'blue')
    title('Normalized Ratio Theta/Delta Power')
    xlabel('Time (s)')
    yline(LIAthresh, 'red')
    yline(runthresh, 'green')
    hold off
    % SLM:
    nexttile % Plot normalized ratio power
    plot(time_spec, ratiosc, 'Color' ,'#dad0d8')
    hold on
    plot(time_spec, ratiofilt, 'blue')
    title('Normalized Ratio Theta/Delta Power')
    xlabel('Time (s)')
    yline(LIAthresh, 'red')
    yline(runthresh, 'green')
    hold off


    % Choose which channel to keep for theta event calculation:


    thetachcalc =  questdlg('Choose channel to compute theta: ', ...
    	'Run-theta events', ...
    	'SP channel','SLM channel','SLM channel');
    switch thetachcalc
        case 'SP channel'
            thetachcalc = 1;
        case 'SLM channel'
            thetachcalc = 2;
    end

    if thetachcalc == 1
        ratiofilt = sp_ratiofilt;
    end

    % Get number of events above that threshold:

    thetabouts = ratiofilt > runthresh;
    [~, locs] = findpeaks(double(thetabouts)); % Para calcular el núm de eventos de theta (realmente te coge el inicio de cada uno de ellos).

   if isempty(locs) %*
    
       results.thetaevent.pscompleteprofile = thetach_profile;
       results.thetaevent.psdelta = thetach_psde;
       results.thetaevent.pstheta = thetach_psth;
       results.thetaevent.thetapstime = nan;
       results.thetaevent.thetapsmean = nan;
       results.thetaevent.gammapstime = nan;
       results.thetaevent.gammapsmean = nan;
       results.thetaevent.eventnr = nan;
       results.thetaevent.evdurtot = nan;
       results.thetaevent.evdurmean = nan;
       results.thetaevent.idxs = nan;
       results.thetaevent.gammahighpsmean = nan;
       results.thetaevent.gammalowpsmean = nan;
       results.thetaevent.spthetafilt = sp_thetafilt;
       results.thetaevent.slmthetafilt = thetafilt;

   
       results.thetaevent.params.runthresh = runthresh; 
       results.thetaevent.params.LIAthresh = LIAthresh;
       if thetachcalc == 1
           results.thetaevent.params.eventch = ripplech;
       elseif thetachcalc == 2
           results.thetaevent.params.eventch = thetach;
       end

       resultsname = char(strcat(lfp_paths(nn), '\LFPdsdata.mat'));
       save(resultsname, 'results');

 

   %% *5. Save data to csv or excel to follow up with stats:
   
   id = results.id.mouseid;
   expcond = results.id.expcond;
   prepost = results.id.time;
   path = results.id.filename;
   eventnr = nan;
   evdurtot = nan;
   evdurmean = nan;
   totmean_theta(nn) = nan;
   totmean_gamma(nn) = nan;
   totmean_gammahigh(nn) = nan;
   totmean_gammalow(nn) = nan;

   Data(nn, :) = {id, expcond, prepost, eventnr, evdurtot, evdurmean, totmean_theta(nn), totmean_gamma(nn),totmean_gammahigh(nn), totmean_gammalow(nn),  path};

   clear events evdur eventnr bouts boutdur betwbouts onset offset shortest subst bouts2
   close all

   else
          % 4.1.1. Get onset and offset of each event saved in locs:

          for kk = 1:length(locs)
              if locs(kk) == locs(end)
                  onset(kk) = locs(kk);
                  if thetabouts(end) == 0
                      offset(kk) = locs(kk) + find(thetabouts(locs(kk):end) == 0, 1, 'first') - 1 ;
                  else
                      offset(kk) = length(thetabouts);
                  end
              else
                 onset(kk) = locs(kk);
                 offset(kk) = locs(kk) + find(thetabouts(locs(kk):end) == 0, 1, 'first') - 1 ;
              end
          end

       % Event matrix (first col is onset, second col is offset):

       bouts = [onset.' offset.'];

       % 4.2. Check if any event is less than 1 second apart from each other, if so,
       % consider it part of the same event (criteria like in Cell 2020)

       betwbouts = bouts(2:end, 1) - bouts(1:end-1, 2);
       subst = betwbouts;
       subst = [subst; 999];
       subst = subst < 2;
       betwbouts = [999; betwbouts];
       shortest = betwbouts < 2;

       if sum(shortest) > 0
           bouts(subst, 2)  = bouts(shortest, 2);
           bouts2(:, :) = bouts(~shortest, :); % Get offset of event closer than <1 sec and assign it to previous event (we interpret it as the same event)
           bouts = bouts2;
       end

       % 4.3. Discard events that are shorter than 1 second:

       boutdur = bouts(:, 2) - bouts(:, 1);
       idx = boutdur < psfs; % non events, get idx (shorter than 1 second)
       if sum(idx) > 0
           bouts = bouts(~idx, :);
       end



    %% 5. GET POWER EVENT MEAN AND TOTAL MEAN OF THETA AND GAMMA
    
    events = bouts;
    evdurtot = sum(boutdur./psfs);
    evdurmean = mean(boutdur./psfs);
    eventnr = size(bouts, 1);
    gammafq = (freqs >= 30 & freqs <= 90); 
    gammalow = (freqs >= 30 & freqs < 60);
    gammahigh = (freqs >= 60 & freqs < 90);

    for jj = 1:eventnr
        varname = sprintf('event%i', jj);
        % theta
        thetaevent_ps(nn).(varname) = thetach_profile(thetafq, events(jj, 1):events(jj, 2));
        thetaevent_pseacheventmean(nn).(varname) = pow2db(mean(mean(thetaevent_ps(nn).(varname), 1), 2));
        th_eachevent(jj) = thetaevent_pseacheventmean(nn).(varname);
        % gamma (all)
        gammaevent_ps(nn).(varname) = thetach_profile(gammafq, events(jj, 1):events(jj, 2)); 
        gammaevent_pseacheventmean(nn).(varname) = pow2db(mean(mean(gammaevent_ps(nn).(varname), 1), 2));
        ga_eachevent(jj) = gammaevent_pseacheventmean(nn).(varname);       
        % gamma(high)
        gammahighevent_ps(nn).(varname) = thetach_profile(gammahigh, events(jj, 1):events(jj, 2)); 
        gammahighevent_pseacheventmean(nn).(varname) = pow2db(mean(mean(gammahighevent_ps(nn).(varname), 1), 2));
        gahigh_eachevent(jj) = gammahighevent_pseacheventmean(nn).(varname);    
        % gamma(low)
        gammalowevent_ps(nn).(varname) = thetach_profile(gammalow, events(jj, 1):events(jj, 2)); 
        gammalowevent_pseacheventmean(nn).(varname) = pow2db(mean(mean(gammalowevent_ps(nn).(varname), 1), 2));
        galow_eachevent(jj) = gammalowevent_pseacheventmean(nn).(varname);    
    end

    totmean_theta(nn) = mean(th_eachevent);
    totmean_gamma(nn) = mean(ga_eachevent);
    totmean_gammahigh(nn) = mean(gahigh_eachevent);
    totmean_gammalow(nn) = mean(galow_eachevent);



   %% 6. Save theta and gamma power (mean and each event separately): 

   results.thetaevent.pscompleteprofile = thetach_profile;
   results.thetaevent.psdelta = thetach_psde;
   results.thetaevent.pstheta = thetach_psth;
   results.thetaevent.thetapstime = thetaevent_ps;
   results.thetaevent.thetapsmean = thetaevent_pseacheventmean;
   results.thetaevent.gammapstime = gammaevent_ps;
   results.thetaevent.gammapsmean = gammaevent_pseacheventmean;
   results.thetaevent.gammahighpsmean = gammahighevent_pseacheventmean;
   results.thetaevent.gammalowpsmean = gammalowevent_pseacheventmean;
   results.thetaevent.spthetafilt = sp_thetafilt;
   results.thetaevent.slmthetafilt = thetafilt;
   
   results.thetaevent.eventnr = eventnr;
   results.thetaevent.evdurtot = evdurtot;
   results.thetaevent.evdurmean = evdurmean;
   results.thetaevent.idxs = events;

   results.thetaevent.params.runthresh = runthresh; % threshold of zscore in std
   results.thetaevent.params.LIAthresh = LIAthresh;
       if thetachcalc == 1
           results.thetaevent.params.eventch = ripplech;
       elseif thetachcalc == 2
           results.thetaevent.params.eventch = thetach;
       end

   resultsname = char(strcat(lfp_paths(nn), '\LFPdsdata.mat'));
   fig1name = char(strcat(lfp_paths(nn), '\SPvSLMps_thpow_ratiopow.jpg'));
   save(resultsname, 'results');
   saveas(fig1, fig1name)


   %% 5. Save data to csv or excel to follow up with stats:

   id = results.id.mouseid;
   expcond = results.id.expcond;
   prepost = results.id.time;
   path = results.id.filename;

   Data(nn, :) = {id, expcond, prepost, eventnr, evdurtot, evdurmean, totmean_theta(nn), totmean_gamma(nn), totmean_gammahigh(nn), totmean_gammalow(nn), path};

   clear events evdur eventnr bouts boutdur betwbouts onset offset shortest subst bouts2
   close all
    end
 end


   %% 6. Save csv files: 
    
       Data = cell2table(Data, 'VariableNames', measurenames);

       savename = uigetdir();
       dataname = strcat(savename, '\tg_thetaevent.xlsx');
       
       writetable(Data, dataname, 'WriteVariableNames', true);


end



      