
function [] = theta_event()


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
   Data1 = Data;



for n = 1:length(lfp_paths)


    %% 1. Load saved file:

    filename = char(strcat(lfp_paths(n), '\LFPdsdata.mat'));
    message = sprintf(['Processing file %i out of %i' ...
        ' ' ...
        ], n, length(lfp_paths));
    disp(message)
    load(filename, 'results')


    %% 2. Extract and generate main variables:

% 2.1. SLM channel:

    thetach = results.params.thetach;
    thetach_profile = results.psProfile.pS(:, :, thetach);
    freqs = results.psProfile.pSfreqs;
    
    %theta_cf = results.psProfile.params.theta_cf;
    theta_cf = [5.5 8.5]; % Mismo rango que Alberto
    thetafq = (results.psProfile.pSfreqs >= theta_cf(1) & results.psProfile.pSfreqs <= theta_cf(2)); % fq range
    thetarg = freqs(thetafq);
    thetach_psth = thetach_profile(thetafq, :); % Espectro sólo para las fq de interés
    
    %delta_cf = results.psProfile.params.delta_cf;
    delta_cf = [1 3.5]; % No puede ser de menos por la ventana de desplazamiento que hemos puesto en power spectrum (2s)
    deltafq = (results.psProfile.pSfreqs >= delta_cf(1) & results.psProfile.pSfreqs <= delta_cf(2)); % based on Malezieux, M., Kees, A. L., & Mulle, C. (2020)
    deltarg = freqs(deltafq);
    thetach_psde= thetach_profile(deltafq, :);

    if size(thetach_psde, 1) > size(thetach_psth, 1)
        thetach_psde = thetach_psde(1:size(thetach_psth, 1), :);
    elseif size(thetach_psde, 1) < size(thetach_psth, 1)
        thetach_psth = thetach_psth(1:size(thetach_psde, 1), :);
    end
    
    psfs = results.psProfile.params.tWin/1000;
%     time_spec = [0:size(thetach_psth, 2)-1]./psfs;
    time_spec = (0:size(thetach_psth, 2)-1)./2; % 2 muestras por seg

% 2.2. SP channel: 

    ripplech = results.params.ripplech;
    ripplech_profile = results.psProfile.pS(:, :, ripplech);
    freqs = results.psProfile.pSfreqs;
    
    %theta_cf = results.psProfile.params.theta_cf;
    sp_thetafq = (results.psProfile.pSfreqs >= theta_cf(1) & results.psProfile.pSfreqs <= theta_cf(2)); % fq range
    ripplech_psth = ripplech_profile(sp_thetafq, :); % Espectro sólo para las fq de interés
%     dbth = pow2db(ripplech_psth);
%     dbde = pow2db(ripplech_psde);
%     dbth2 = pow2db(thetach_psth);
%     dbde2 = pow2db(thetach_psde);
    
    %delta_cf = results.psProfile.params.delta_cf;
    sp_deltafq = (results.psProfile.pSfreqs >= delta_cf(1) & results.psProfile.pSfreqs <= delta_cf(2)); % based on Malezieux, M., Kees, A. L., & Mulle, C. (2020)
    ripplech_psde= ripplech_profile(sp_deltafq, :);

    if size(ripplech_psde, 1) > size(ripplech_psth, 1)
        ripplech_psde = ripplech_psde(1:size(ripplech_psth, 1), :);
    elseif size(ripplech_psde, 1) < size(ripplech_psth, 1)
        ripplech_psth = ripplech_psth(1:size(ripplech_psde, 1), :);
    end
    
    psfs = 2; % Pendiente de automatizar... no sé por qué no está saliendo lo que debería
%     psfs = results.psProfile.params.tWin/1000;
%     time_spec = [0:size(thetach_psth, 2)-1]./psfs;
    time_spec = (0:size(ripplech_psth, 2)-1)./2; % 2 muestras por seg


    
    %% 3. Ratio theta/delta:

% 3.1. SLM channel:
thetamean = mean(thetach_psth, 1); % Vector de tiempo >> media de fq
deltamean = mean(thetach_psde, 1);
ratiotd = thetamean./deltamean; 
% ratiotd = deltamean./thetamean;
if sum(isinf(ratiotd)) > 0
    ratiotd(ratiotd == Inf) = nan;
end

% Código de Alberto para sgolayfilt:
ratiofilt = movmean( movmean( sgolayfilt(ratiotd,50,51), 20), 20); % Sale igual que si lo normalizas antes
ratiofilt = rescale(ratiofilt, 0, 1); 
% ratiofiltmean = mean(ratiofilt, 'omitnan');
% ratiofiltstd = std(ratiofilt);

ratiosc = rescale(ratiotd, 0, 1); % normalización como Alberto a 0-1. 
% ratioscmean = mean(ratiosc, 'omitnan');
% ratioscstd = std(ratiosc, 'omitnan');

% Apply sgolayfilt + rescale also to theta power (como Alberto): 

thetameansc = rescale(thetamean, 0, 1);
thetafilt =  movmean( movmean( sgolayfilt(thetamean,50,51), 20), 20);
thetafilt = rescale(thetafilt, 0, 1);


% 3.2. SP channel: 

sp_thetamean = mean(ripplech_psth, 1); % Vector de tiempo >> media de fq
sp_deltamean = mean(ripplech_psde, 1);
sp_ratiotd = sp_thetamean./sp_deltamean; 
if sum(isinf(sp_ratiotd)) > 0
    sp_ratiotd(sp_ratiotd == Inf) = nan;
end

sp_ratiosc = rescale(sp_ratiotd, 0, 1); % normalización como Alberto a 0-1. 
% sp_ratioscmean = mean(sp_ratiosc, 'omitnan');
% sp_ratioscstd = std(sp_ratiosc, 'omitnan');

% Código de Alberto para sgolayfilt:
sp_ratiofilt = movmean( movmean( sgolayfilt(sp_ratiotd,50,51), 20), 20);
sp_ratiofilt = rescale(sp_ratiofilt, 0, 1); 
sp_ratiofiltmean = mean(sp_ratiofilt, 'omitnan');
sp_ratiofiltstd = std(sp_ratiofilt);

sp_thetameansc = rescale(sp_thetamean, 0, 1);
sp_thetafilt =  movmean( movmean( sgolayfilt(sp_thetamean,50,51), 20), 20);
sp_thetafilt = rescale(sp_thetafilt, 0, 1);




%%%%%%%%%%%%%%%% CRITERIA TO DETECT RUN-THETA EVENTS %%%%%%%%%%%%%%%%%%%%%
% 1. Thresholds: 0.5 for run theta. This threshold is applied to filtered
% ratio (env, sgolayfilt). 
% 2. Events that are less than 1 sec apart from each other are considered
% part of the same event. 
% 3. Events shorter than 1 sec are not considered as theta events

% 1.  THRESHOLDS COMO ALBERTO:

LIAthresh = 0.2; % save to detect ripples events under this threshold
runthresh = 0.3; % Estos valores después de normalizar entre 0 y 1

% Compare SP and SLM to check which channel to choose for theta calc.: 

    imrg = results.psProfile.pSfreqs >= 1 & results.psProfile.pSfreqs <= 15;
    imgfreqs = freqs(imrg);

fig1 = figure(n);

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
        thetachcalc = 0;
end

if thetachcalc == 0
    ratiofilt = sp_ratiofilt;
end

% Get number of events above that threshold:

thetabouts = ratiofilt > runthresh;
[~, locs] = findpeaks(double(thetabouts)); % Para calcular el núm de eventos de theta (realmente te coge el inciio de cada uno de ellos). 

   if isempty(locs)
    
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
   
       results.thetaevent.params.runthresh = runthresh; % threshold of zscore in std
       results.thetaevent.params.LIAthresh = LIAthresh;
       if thetachcalc == 1
           results.thetaevent.params.eventch = ripplech;
       elseif thetachcal == 0
           results.thetaevent.params.eventch = thetach;
       end

       resultsname = char(strcat(lfp_paths(n), '\LFPdsdata.mat'));
%        figname1 = char(strcat(lfp_paths(n), '\psthetadelta.jpg'));
%        figname2 = char(strcat(lfp_paths(n), '\ZscoreRatioTD.jpg'));
       save(resultsname, 'results');
%        saveas(fig1, figname1)
%        saveas(fig2, figname2)


   %% 5. Save data to csv or excel to follow up with stats:
   
   id = results.id.mouseid;
   expcond = results.id.expcond;
   prepost = results.id.time;
   path = results.id.filename;
   eventnr = nan;
   evdurtot = nan;
   evdurmean = nan;
   totmean_theta(n) = nan;
   totmean_gamma(n) = nan;
    totmean_gammahigh(n) = nan;
    totmean_gammalow(n) = nan;

   Data1(n, :) = {id, expcond, prepost, eventnr, evdurtot, evdurmean, totmean_theta(n), totmean_gamma(n),totmean_gammahigh(n), totmean_gammalow(n),  path};

   clear events evdur eventnr bouts boutdur betwbouts onset offset shortest subst bouts2
   close all

   else
       % Get onset and offset of each event saved in locs:

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

       % 2. Check if any event is less than 1 second apart from each other, if so,
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

       % 3. Discard events that are shorter than 1 second:

       boutdur = bouts(:, 2) - bouts(:, 1);
       idx = boutdur < psfs; % non events, get idx (shorter than 1 second)
       if sum(idx) > 0
           bouts = bouts(~idx, :);
       end



    %% GET POWER EVENT MEAN AND TOTAL MEAN OF THETA AND GAMMA
    
    events = bouts;
    evdurtot = sum(boutdur./psfs);
    evdurmean = mean(boutdur./psfs);
    eventnr = size(bouts, 1);
    gammafq = (freqs >= 30 & freqs <= 90); % Ojo que para esto habría
    % que aplicar filtro antes a los datos debido a la señal de ruido
    % eléctrico que se genera en 50 Hz 
    gammalow = (freqs >= 30 & freqs < 60);
    gammahigh = (freqs >= 60 & freqs < 90);

%  if thetachcalc == 1 % SP channel
%     for jj = 1:eventnr
%         varname = sprintf('event%i', jj);
%         % theta
%         thetaevent_ps(n).(varname) = ripplech_profile(thetafq, events(jj, 1):events(jj, 2));
%         thetaevent_pseacheventmean(n).(varname) = pow2db(mean(mean(thetaevent_ps(n).(varname), 1), 2));
%         th_eachevent(jj) = thetaevent_pseacheventmean(n).(varname);
%         % gamma (all)
%         gammaevent_ps(n).(varname) = ripplech_profile(gammafq, events(jj, 1):events(jj, 2)); 
%         gammaevent_pseacheventmean(n).(varname) = pow2db(mean(mean(gammaevent_ps(n).(varname), 1), 2));
%         ga_eachevent(jj) = gammaevent_pseacheventmean(n).(varname);       
%         % gamma(high)
%         gammahighevent_ps(n).(varname) = ripplech_profile(gammahigh, events(jj, 1):events(jj, 2)); 
%         gammahighevent_pseacheventmean(n).(varname) = pow2db(mean(mean(gammahighevent_ps(n).(varname), 1), 2));
%         gahigh_eachevent(jj) = gammahighevent_pseacheventmean(n).(varname);    
%         % gamma(low)
%         gammalowevent_ps(n).(varname) = ripplech_profile(gammalow, events(jj, 1):events(jj, 2)); 
%         gammalowevent_pseacheventmean(n).(varname) = pow2db(mean(mean(gammalowevent_ps(n).(varname), 1), 2));
%         galow_eachevent(jj) = gammalowevent_pseacheventmean(n).(varname);    
%     end

% elseif thetachcalc == 0 % SLM channel
    for jj = 1:eventnr
        varname = sprintf('event%i', jj);
        % theta
        thetaevent_ps(n).(varname) = thetach_profile(thetafq, events(jj, 1):events(jj, 2));
        thetaevent_pseacheventmean(n).(varname) = pow2db(mean(mean(thetaevent_ps(n).(varname), 1), 2));
        th_eachevent(jj) = thetaevent_pseacheventmean(n).(varname);
        % gamma (all)
        gammaevent_ps(n).(varname) = thetach_profile(gammafq, events(jj, 1):events(jj, 2)); 
        gammaevent_pseacheventmean(n).(varname) = pow2db(mean(mean(gammaevent_ps(n).(varname), 1), 2));
        ga_eachevent(jj) = gammaevent_pseacheventmean(n).(varname);       
        % gamma(high)
        gammahighevent_ps(n).(varname) = thetach_profile(gammahigh, events(jj, 1):events(jj, 2)); 
        gammahighevent_pseacheventmean(n).(varname) = pow2db(mean(mean(gammahighevent_ps(n).(varname), 1), 2));
        gahigh_eachevent(jj) = gammahighevent_pseacheventmean(n).(varname);    
        % gamma(low)
        gammalowevent_ps(n).(varname) = thetach_profile(gammalow, events(jj, 1):events(jj, 2)); 
        gammalowevent_pseacheventmean(n).(varname) = pow2db(mean(mean(gammalowevent_ps(n).(varname), 1), 2));
        galow_eachevent(jj) = gammalowevent_pseacheventmean(n).(varname);    
    end
% end
    totmean_theta(n) = mean(th_eachevent);
    totmean_gamma(n) = mean(ga_eachevent);
    totmean_gammahigh(n) = mean(gahigh_eachevent);
    totmean_gammalow(n) = mean(galow_eachevent);



   %% 4. Save theta and gamma power (mean and each event separately): 

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
       elseif thetachcalc == 0
           results.thetaevent.params.eventch = thetach;
       end

   resultsname = char(strcat(lfp_paths(n), '\LFPdsdata.mat'));
   fig1name = char(strcat(lfp_paths(n), '\SPvSLMps_thpow_ratiopow.jpg'));
%    figname1 = char(strcat(lfp_paths(n), '\psthetadelta.jpg'));
%    figname2 = char(strcat(lfp_paths(n), '\ZscoreRatioTD.jpg'));
   save(resultsname, 'results');
   saveas(fig1, fig1name)
%    saveas(fig1, figname1)
%    saveas(fig2, figname2)


   %% 5. Save data to csv or excel to follow up with stats:

   id = results.id.mouseid;
   expcond = results.id.expcond;
   prepost = results.id.time;
   path = results.id.filename;

   Data1(n, :) = {id, expcond, prepost, eventnr, evdurtot, evdurmean, totmean_theta(n), totmean_gamma(n), totmean_gammahigh(n), totmean_gammalow(n), path};

   clear events evdur eventnr bouts boutdur betwbouts onset offset shortest subst bouts2
   close all
       end
 end


   %% 6. Save csv files: 
    
       Data1 = cell2table(Data1, 'VariableNames', measurenames);

       savename = uigetdir();
       data1name = strcat(savename, '\tg_thetaevent.xlsx');
       
       writetable(Data1, data1name, 'WriteVariableNames', true);


end



      