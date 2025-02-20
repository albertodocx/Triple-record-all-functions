
function [] = getripples()

%% 0. Choose folder and load filepaths

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

% Prepare variables for data analysis:
id = '';
expcond = '';
prepost = '';
thetaprop = [];
thetaevpermin = [];
LIAprop = [];
LIAevpermin = [];
meanrippleqt = [];
totalrippleqt = [];
meanLIAdur = [];
totalLIAdur = [];
meanripplefq_table = [];
meanripdur_table = [];
meansdpowrip_table = [];
meanmaxpowrip_table = [];
path = '';

measurenames = {'Id' 'expcond' 'prepost' 'Theta Time(%)' 'Theta Event per min(n/min)' 'LIA Time (%)' 'LIA Event per min (n/min)' 'mean rippleqt (n/LIA event)' 'total ripple qt (n)' 'mean LIA event duration (s)' 'total LIA event duration (s)' 'Mean Ripple Fq (Hz)' 'Mean Ripple Duration (ms)' 'Mean Power (SD)' 'Max Power (SD)' 'path'};
Data = {id, expcond, prepost,thetaprop, thetaevpermin, LIAprop, LIAevpermin, meanrippleqt, totalrippleqt, meanLIAdur, totalLIAdur, meanripplefq_table ,meanripdur_table, meansdpowrip_table,meanmaxpowrip_table,  path};
Data1 = Data;


for nn = 1:length(lfp_paths)


    %% 1. Load saved file and other variables:

    filename = char(strcat(lfp_paths(nn), '\LFPdsdata.mat'));
    message = sprintf(['Processing file %i out of %i' ...
        ' ' ...
        ], nn, length(lfp_paths));
    disp(message)
    load(filename, 'results')

    % 1.1. Load variables:

    % Channels and ripple channel data:
    ripplech = results.params.ripplech;
    thetach = results.params.thetach;
    rippledata = results.LFPds.dsdata(:, ripplech);
    LFPdata = results.LFPds.dsdata;
    
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
    LIAthresh = results.thetaevent.params.LIAthresh;

    % Theta/delta ratios for each channel: 
    spratio = results.thetaevent.spthetafilt;
    slmratio = results.thetaevent.slmthetafilt;

    % Generate time variable for LFP: 
    
    ts = (0:length(rippledata)-1)./lfpfs;
    psts = (0:length(sprippleps)-1)./psfs;
    
    % Choose same ratio as the one used for theta event:
    
    eventch = results.thetaevent.params.eventch;
    

    %% 2. Get ripple events: 
    
    if eventch == ripplech
        LIAidx = spratio < LIAthresh;
    elseif eventch == thetach
        LIAidx = spratio < LIAthresh;
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
    
    fctr = lfpfs/psfs;

    if LIAbouts(1, 1) == 1
        LIAbouts(2:end, 1) = LIAbouts(2:end, 1)*fctr;
        LIAbouts(:, 2) = LIAbouts(:, 2) * fctr;
    else
        LIAbouts = LIAbouts * fctr;
    end
%prueba
%     tiledlayout(2, 4)
% for ii = 1:length(locs)
%     nexttile
%     plot(ts(LIAbouts(ii, 1):LIAbouts(ii, 2)), rippledata(LIAbouts(ii, 1):LIAbouts(ii, 2)));
% end


    %% 3. Apply bandpass filter to raw signal - usar filtroFIR1.m (de Alberto): 

    %     for ii = 1:length(locs)
    %         rippeventsdata{ii, 1} = rippledata(ripplebouts(ii, 1):ripplebouts(ii, 2));
    %         rippeventsdata{ii, 2} = bandpass(rippeventsdata{ii, 1}, [90 200], lfpfs);
    %         rippeventsdata{ii, 3} =  movmean( movmean( sgolayfilt(rippeventsdata{ii, 2},50,51), 20), 20);
    %     end

    % Of the whole channel and LFPdata: 

    pass1 = 90; % Hz
    pass2 = 200; % Hz
    
    [signRip_complete, signEnv_complete] = filtroFIR1(lfpfs, LFPdata, pass1, pass2);
    [signRip_ripples, signEnv_ripples] = filtroFIR1(lfpfs, rippledata, pass1, pass2);
   
    ripplethresh = mean(signEnv_ripples + 3*std(signEnv_ripples));

%     %% PRUEBAS PARA COMPROBAR CÓMO SE COMPORTA EL FILTRO
%     plot(ts(LIAbouts(1):LIAbouts(2)), signRip_ripples(LIAbouts(1):LIAbouts(2)),  'Color' ,'#dad0d8')
%     hold on
%     yline(ripplethresh)
%     plot(ts(LIAbouts(1):LIAbouts(2)), signEnv_ripples(LIAbouts(1):LIAbouts(2)), 'blue')
%     hold off
% 
%     tiledlayout(2, 1)
%     nexttile
%     plot(ts(LIAbouts(1):LIAbouts(2)), signEnv_ripples(LIAbouts(1):LIAbouts(2)), 'blue')
%     hold on
%     yline(ripplethresh)
%     hold off
%     nexttile
%     plot(ts(LIAbouts(1):LIAbouts(2)), double(tempidx))
%     
%      plot(rippledur.durinms)
%      hold on
%      yline(150)
%      yline(40)
% 

    %% 4. Isolate ripple bouts:

    ripplebouts = struct();

   for ii = 1:length(locs)

    temprip = signEnv_ripples(LIAbouts(ii, 1):LIAbouts(ii, 2)); % envelope of each LIA event
    tempts = ts(LIAbouts(ii, 1):LIAbouts(ii, 2));
    tempidx = temprip > ripplethresh; % Find idx above that threshold (full vector)
    % How many ripples and where are they?
    [~, ripidx] = findpeaks(double(tempidx));        
    if tempidx(1) == 1
        ripidx = [1; ripidx]; % ripple event
    end
    
           % Get onset and offset of each event saved in locs:
       if ~isempty(ripidx)
           for kk = 1:length(ripidx)
               if ripidx(kk) == ripidx(end)
                   temponset(kk) = ripidx(kk);
                   if tempidx(end) == 0
                       tempoffset(kk) = ripidx(kk) + find(tempidx(ripidx(kk):end) == 0, 1, 'first') - 1 ;
                   else
                       tempoffset(kk) = length(tempidx);
                   end
               else
                   temponset(kk) = ripidx(kk);
                   tempoffset(kk) = ripidx(kk) + find(tempidx(ripidx(kk):end) == 0, 1, 'first') - 1 ;
               end
           end
           ripplebouts(ii).event = [temponset.', tempoffset.'];
           rippleqt(ii) = length(ripidx);
           rippledur(ii).durinsec = (tempoffset - temponset)./lfpfs;
           rippledur(ii).durinms = rippledur(ii).durinsec * 1000;
           LIAeventdur(ii) = length(temprip)./lfpfs; % En segundos
           envrip(ii).event = temprip;
           filtrip(ii).event = signRip_ripples(LIAbouts(ii, 1):LIAbouts(ii, 2));
           lfpcheck(ii).event = LFPdata(LIAbouts(ii, 1):LIAbouts(ii, 2), :);
           tsrip(ii).event = tempts;

       elseif isempty(ripidx)
           ripplebouts(ii).event = nan;
           rippleqt(ii) = 0;
           LIAeventdur(ii) = length(temprip)./lfpfs; % En segundos      
           envrip(ii).event = temprip;
           filtrip(ii).event = signRip_ripples(LIAbouts(ii, 1):LIAbouts(ii, 2));
           lfpcheck(ii).event = LFPdata(LIAbouts(ii, 1):LIAbouts(ii, 2), :);
           tsrip(ii).event = tempts;               
       end
    clear temprip tempidx ripidx temponset tempoffset 

   end

   %% 5. Select manually which ripples to keep: 

   % Plot each ripple to choose which ones we keep:

   jump = [2, 4, 6, 8, 10, 12, 14, 16];
   envmean = mean(signEnv_ripples);
   envsd = std(signEnv_ripples);

   for ii = 1:length(locs)
       for jj = 1:rippleqt(ii)
           if rippleqt(ii) == 0
               continue
           else
               if (ripplebouts(ii).event(jj,1)-0.2*lfpfs) < 0
                   tempts = tsrip(ii).event(1:ripplebouts(ii).event(jj, 2)+0.2*lfpfs); % menor que longitud de LIA event
               elseif (ripplebouts(ii).event(jj, 2)+0.2*lfpfs) > length(tsrip(ii).event)
                   tempts = tsrip(ii).event(ripplebouts(ii).event(jj,1)-0.2*lfpfs:ripplebouts(ii).event(jj, 2)); % menor q longitud de LIA event                   
               else
                   tempts = tsrip(ii).event(ripplebouts(ii).event(jj,1)-0.2*lfpfs:ripplebouts(ii).event(jj, 2)+0.2*lfpfs);
               end
               tempenvrip = envrip(ii).event(ripplebouts(ii).event(jj,1):ripplebouts(ii).event(jj, 2)); % sólo el ripple


               if rippledur(ii).durinms(jj) < 30 || rippledur(ii).durinms(jj) > 150 % Descartados por duración 
                       ripplestokeep.rippledur(ii).durinms(jj) = nan;
                       ripplestokeep.meansd(ii).event(jj) = nan; % mean sd of all sd values above threshold
                       ripplestokeep.maxsd(ii).event(jj) = nan; % sd of max of envelope
                       continue
               else
               figure('Name','Keep Ripple?', 'Position',[10 10 600 900] ) 
                   for kk = 1:8
                       if (ripplebouts(ii).event(jj,1)-0.2*lfpfs) < 0
                           templfp = lfpcheck(ii).event(1:ripplebouts(ii).event(jj, 2)+0.2*lfpfs, jump(kk));
                       elseif  (ripplebouts(ii).event(jj, 2)+0.2*lfpfs) > length(tsrip(ii).event)
                           templfp = lfpcheck(ii).event(ripplebouts(ii).event(jj,1)-0.2*lfpfs:ripplebouts(ii).event(jj, 2), jump(kk));
                       else
                           templfp = lfpcheck(ii).event(ripplebouts(ii).event(jj,1)-0.2*lfpfs:ripplebouts(ii).event(jj, 2)+0.2*lfpfs, jump(kk));
                       end
                       subplot(8, 1, kk)
                       plot(tempts, templfp)
                       hold on
                       x_points = [tsrip(ii).event(ripplebouts(ii).event(jj,1)), tsrip(ii).event(ripplebouts(ii).event(jj,1)), tsrip(ii).event(ripplebouts(ii).event(jj,2)), tsrip(ii).event(ripplebouts(ii).event(jj,2))];
                       y_points = [min(templfp), max(templfp), max(templfp), min(templfp)];
                       color = [0, 0, 1];
                       a = fill(x_points, y_points, color);
                       a.FaceAlpha = 0.1;
                       channelname = sprintf('channel%i', jump(kk));
                       title(channelname)
                       hold off
                   end

                   questmessage = sprintf('Keep ripple? (Pyrchannel is %i)', ripplech);
                   keepripples = questdlg(questmessage, ...
                     	'Ripples Check', ...
                    	'Yes','No','No');
                   % Handle response
                   switch keepripples
                       case 'Yes'
                           % Sacar medidas por evento (después se hace la media, pero
                           % la matriz es más facil de manipular que las tablas
                           % dinámicas)
                           if rippledur(ii).durinms(jj) >= 30 || rippledur(ii).durinms(jj) <= 150
                               ripplestokeep.rippledur(ii).durinms(jj) = rippledur(ii).durinms(jj);
                               ripplestokeep.meansd(ii).event(jj) = mean((tempenvrip-envmean)/envsd); % mean sd of all sd values above threshold
                               ripplestokeep.maxsd(ii).event(jj) = (max(tempenvrip) - envmean)/envsd; % sd of max of envelope
                           end
                       case 'No'
                           ripplestokeep.rippledur(ii).durinms(jj) = nan;
                           ripplestokeep.meansd(ii).event(jj) = nan; % mean sd of all sd values above threshold
                           ripplestokeep.maxsd(ii).event(jj) = nan; % sd of max of envelope
                   end
               end
           end
       end
           % Save LIA events, with other vars in the same struct >>> not necessary
           %     ripplestokeep.envrip(ii).event(jj) = envrip(ii).event(jj);
           %     ripplestokeep.filtrip(ii).event(jj) = filtrip(ii).event(jj);
           %     ripplestokeep.tsrip(ii).event(jj) = tsrip(ii).event(jj);
           % Aquí calcular la fq media del ripple durante un evento de LIA:

           if rippleqt(ii) == 0
               ripplestokeep.rippleqt(ii) = nan;
               meanripplefq(ii) = nan;           
           else 
               ripplestokeep.rippleqt(ii) = sum(~isnan(ripplestokeep.rippledur(ii).durinms));
               meanripplefq(ii) = ripplestokeep.rippleqt(ii)./LIAeventdur(ii);
           end
       close all
   end

   %% 6. Calculate and arrange ripple measurements:

   % 6.1. Each LIA event:

   for ii = 1:length(locs)-1
       meanrippledur(ii) = mean(ripplestokeep.rippledur(ii).durinms, 'omitnan'); % Por evento
       meanripplesd(ii) = mean(ripplestokeep.meansd(ii).event, 'omitnan');
       maxripplesd(ii) = mean(ripplestokeep.maxsd(ii).event, 'omitnan');
   end

   % 6.2. Reduce all values to one per recording (mean of means): 

       totmeanrippledur = mean(meanrippledur, 'omitnan');
       totmeanripplesd = mean(meanripplesd, 'omitnan');
       totmaxripplesd = mean(maxripplesd, 'omitnan');


       %% 7. Save data:

   id = results.id.mouseid;
   expcond = results.id.expcond;
   prepost = results.id.time;
   path = results.id.filename;
   meanrippleqt = mean(ripplestokeep.rippleqt, 'omitnan'); % incluye 0s
   totalrippleqt = sum(ripplestokeep.rippleqt, 'omitnan');
   meanLIAdur = mean(LIAeventdur, 'omitnan');
   totalLIAdur = sum(LIAeventdur, 'omitnan');
   meanripplefq_table = mean(meanripplefq, 'omitnan');
   meanripdur_table = totmeanrippledur;
   meansdpowrip_table = totmeanripplesd;
   meanmaxpowrip_table = totmaxripplesd;
   %% Theta/LIA proportion of time:
   
   thetaprop = results.thetaevent.evdurtot/results.params.durinsec; % Percentage time of theta during recording
   thetaevpermin = results.thetaevent.eventnr/results.params.durinmin; % Events per min of theta during recording
   LIAprop = totalLIAdur/ results.params.durinsec;
   LIAevpermin = size(LIAbouts, 1)/results.params.durinmin ;


   Data1(nn, :) = {id, expcond, prepost,thetaprop, thetaevpermin, LIAprop, LIAevpermin, meanrippleqt, totalrippleqt, meanLIAdur, totalLIAdur, meanripplefq_table ,meanripdur_table, meansdpowrip_table,meanmaxpowrip_table,  path};

   results.rippleevent.ripplestokeep = ripplestokeep;
   results.rippleevent.signEnv_complete = signEnv_complete;
   results.rippleevent.signEnv_ripples = signEnv_ripples;
   results.rippleevent.signRip_complete = signRip_complete;
   results.rippleevent.signRip_ripples = signRip_ripples;
   results.rippleevent.envmean = envmean;
   results.rippleevent.envsd = envsd;
   results.rippleevent.envrip = envrip;
   results.LIAparams.LIAbouts = LIAbouts;
   results.LIAparams.LIAeventdur = LIAeventdur;
   results.LIAparams.ripplethresh = ripplethresh;

   results.thetaevent.thetaprop = thetaprop;
   results.thetaevent.thetaevpermin = thetaevpermin;
   results.rippleevent.LIAprop = LIAprop;
   results.rippleevent.LIAevpermin = LIAevpermin;
   

end % Del bucle por cada registro

   %% 8. Save csv files: 
    
       Data1 = cell2table(Data1, 'VariableNames', measurenames);

       savename = uigetdir();
       data1name = strcat(savename, '\rippledata.xlsx');
       
       writetable(Data1, data1name, 'WriteVariableNames', true);




end % de la función
