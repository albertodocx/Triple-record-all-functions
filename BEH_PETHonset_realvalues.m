function [results] = BEH_PETHonset_realvalues(contvar, fscontvar, evsamp, fsev, combiname,  savepath, varargin)

%% PETHonset plots and analyzes beginning of behavior events

% INPUTS:
%   contvar (vector): continuous variable that you want to allign with specific
%  events (z values). DFF, velocity, acceleration, etc. 
%   fscontvar (escalar): frequency sampling of contvar
%   evsamp (two col matrix): matrix of two columns containing onset sample and offset sample of such
%  events. This could be: freezing, shocks, movement, etc. 
%   fsev (escalar): frequency sampling of events, needed as the fs could be different
%  for those two variables. 
%   savepath (string): where to save all the data

% OPTIONAL: 
%  Pre (s): start of PETH window. Ex: 2 is -2 seconds before event. 
%  Post (s): end of PETH window. 2 is duration of event 2 seconds after
%  start
%  bin: size of each square in the PETH
%  threshold (s): criteria to eliminate events (or not show them)
%  AUCqt (n): how many AUC to calculate
%  AUCint (s): intervals to get information from (coherent with PETH
%  window size)
%  summarymeasures (string): which measures to extract for later
%  quantitative analysis.
%  maxlength (s): max length of your recording/file/etc.  


p = inputParser;
addParameter(p,'Pre', 2,@isnumeric);
addParameter(p,'Post', 2, @isnumeric);
addParameter(p,'bin', 0.1, @isnumeric);
addParameter(p,'threshold', 0.5, @isnumeric);
addParameter(p,'AUCqt', 2, @isnumeric);
addParameter(p,'AUCint', [-2 0; 0 2], @ismatrix);
addParameter(p,'summarymeasures', ["mean", "std"], @isstring);
addParameter(p,'maxlength', 300 ,@isnumeric);


parse(p,varargin{:});
Pre = p.Results.Pre;
Post = p.Results.Post;
bin = p.Results.bin;
threshold = p.Results.threshold;
AUCqt = p.Results.AUCqt;
AUCint = p.Results.AUCint;
summarymeasures = p.Results.summarymeasures;
maxlength = p.Results.maxlength;


%% 1. Load data: 

Events = evsamp./fsev; 
tscontvar = round(Events.*fscontvar);
Time = (0:length(contvar)-1)./fscontvar;

if Events(end, 2) > maxlength % Esto no es generalizable.. 
    Events(end, 2) =  maxlength; % en seg
    tscontvar(end, 2) = length(contvar);
end



%% 2. Variables for PETH: 

 binsize = round(bin*fsev); % Unit is in samples, if each bin represents 0.2 seconds, then how many samples are within 0.2 seconds (in FP fsev)
 PreWind = round(Pre*fsev);
 PostWind = round(Post*fsev);
%  BLWinStart = BinWinStart;
%  BLWinDur = BinWinDuration;
%  binsizeperev = repelem(binsize, size(Events, 1)).';

%% 3. Save which events have overlapping windows

evdif = Events(2:end, 1) - Events(1:end-1, 2); % Diferencia entre eventos (idx +1)
overlap =  evdif < Pre;
% checkfirst = Events(1, 1) - Pre;
% checklast = Events(end, 2) + Post;


%% 4. Check if first event is too close to the beginning and needs also binsize correction

% if checkfirst < Pre
%     binsizeperev(1) = round(Events(1, 1)*fsev); % es esta fsev porque el binsize se coge en función de BHD. 
% end

% OJO QUE BINSIZE ES EL MISMO, LO QUE CAMBIA ES LA LONGITUD DEL BASELINE
% QUE SE CALCULA. 


%% 5. Max length of event >> esto es la PostWind: 

boutdur = Events(:, 2) - Events(:, 1);
[val, ~] = max(boutdur); Postwfr = val; PostWind = round(Postwfr*fsev);

BLdata = zeros([size(Events, 1) PreWind]);
EVdata = zeros([size(Events, 1) PostWind]);


%% 5. Extract FP data from each event and save data:

% COMPLETE DATA WITH OVERLAP

for ii = 1:size(Events, 1)
    if ii == 1
        EVdata(ii, :) = contvar(tscontvar(ii, 1):tscontvar(ii, 1)+ PostWind -1);
        if (tscontvar(ii, 1) - PreWind - 1) < 0 % que exceda por deltante
            % añadir tantos nans por delante == desde inicio de freezing hasta
            % distancia que haya para llegar a PreWind muestras (osea PreWind -
            % las muestras que haya en ese intervalo >> boutdurTS) PreWind -
            % boutdurTS
            BLdata(ii, :) = [nan([1 abs(tscontvar(ii, 1) - PreWind - 1)]) contvar(1:tscontvar(ii, 1)-1)];
        else % caso normal
            BLdata(ii, :) = contvar(tscontvar(ii, 1)-PreWind:tscontvar(ii, 1)-1);
        end
    elseif ii > 1 && ii < size(Events, 1)
        if (tscontvar(ii, 1) - PreWind - 1) < 0 % Que exceda por delante
            BLdata(ii, :) = [nan([1 abs(tscontvar(ii, 1) - PreWind - 1)]) contvar(1:tscontvar(ii, 1)-1)];
            EVdata(ii, :) = contvar(tscontvar(ii, 1):tscontvar(ii, 1)+ PostWind-1);
        elseif (tscontvar(ii, 1) + PostWind) > length(contvar) % Que exceda por detrás
            nansize = (PostWind - length(contvar(tscontvar(ii, 1):end))) ;
            BLdata(ii, :) = contvar(tscontvar(ii, 1)-PreWind:tscontvar(ii, 1)-1);
            EVdata(ii, :) = [contvar(tscontvar(ii, 1):end).', nan([1 nansize ])];
            % añadir nans por el final
        else % caso normal
            BLdata(ii, :) = contvar(tscontvar(ii, 1)-PreWind:tscontvar(ii, 1)-1);
            EVdata(ii, :) = contvar(tscontvar(ii, 1):tscontvar(ii, 1)+ PostWind -1);
        end

    elseif ii == size(Events, 1)
        if (tscontvar(ii, 1) + PostWind -1) > length(contvar)
            nansize = (PostWind - length(contvar(tscontvar(ii, 1):end))) ;
            BLdata(ii, :) = contvar(tscontvar(ii, 1)-PreWind:tscontvar(ii, 1)-1);
            EVdata(ii, :) = [contvar(tscontvar(ii, 1):end).', nan([1 nansize ])];
        else 
            BLdata(ii, :) = contvar(tscontvar(ii, 1)-PreWind:tscontvar(ii, 1)-1);
            EVdata(ii, :) = contvar(tscontvar(ii, 1):tscontvar(ii, 1)+ PostWind -1);
        end
    end
end




%% 6. Correct for length: assign nans where there is no event
% BASELINE: put nans for overlapping trials: (first event is already
% corrected)

overlapsize = zeros([size(overlap) 1]);
overlapsize(overlap) = evdif(overlap) - Pre;
overlapsize = [0; overlapsize];
overlapsize = abs(round(overlapsize.*fsev));

nanBLdata = BLdata;


for jj = 1:size(BLdata, 1)
    if overlapsize(jj) == 0
        continue
    else
    nanvect = nan([1 overlapsize(jj)]);
    nanBLdata(jj, 1:overlapsize(jj)) = nanvect;
    end
end


% EVENT:

nanEVdata = EVdata; % Estas vars ahora contienen nan, de ahí el nombre

boutdursamp = round(boutdur.*fsev);
elimsample = size(EVdata, 2) - boutdursamp;

for kk = 1:size(EVdata, 1) 
    nanvect = nan([1 elimsample(kk)]);
    nanEVdata(kk, boutdursamp(kk)+1:end) = nanvect;
end


%% 7. Compute continuous variable

% for nn = 1:size(BLdata)
%     BLmed(nn) = median(nanBLdata(nn, :), 'omitnan'); BLmad(nn) = mad(nanBLdata(nn, :));
%     EVdataZ(nn, :) = (nanEVdata(nn, :) - BLmed(nn))./BLmad(nn);  % contvar on baseline
% % WITH NORMALIZATION OF BL:
%     BLdataZ(nn, :) = (nanBLdata(nn, :) - BLmed(nn))./BLmad(nn);
% end


%% 8. Downsample through binsize for PETH:

BLjumps = 1:binsize:size(nanBLdata, 2);
EVjumps = 1:binsize:size(nanEVdata, 2);

for ii = 1:length(BLjumps)-1
    binnedBL(:,ii) = median(nanBLdata(:, BLjumps(ii):BLjumps(ii+1)), 2, 'omitnan');
end


for ii = 1:length(EVjumps)-1
    binnedEV(:, ii) = median(nanEVdata(:, EVjumps(ii):EVjumps(ii+1)), 2, 'omitnan');
end

%% 9. Prepare Data for PETH:


% Criteria for BL (threshold = 0.5 seg):

lenBL = zeros([size(nanBLdata, 1), 1]);
for kk = 1:size(nanBLdata, 1)
    lenBL(kk) = sum(~isnan(nanBLdata(kk, :)));
end

sampthresh = round(threshold*size(nanBLdata, 2)./Pre);
underthresh = lenBL < sampthresh;

% Data and plot PETH:

PETHdata = [binnedBL binnedEV];
PETHev = PETHdata; % Save to other var to be saved later wo the criteria and full data

PETHdata = PETHdata(~underthresh, :);
PETHts = linspace(-Pre, Postwfr, size(PETHdata, 2)); 
PETHevnr = size(PETHdata, 1);

% Extract summary measures:

if summarymeasures(1) == "mean"
    PETHm = mean(PETHev, 1, 'omitnan');
elseif summarymeasures(1) == "median"
    PETHm = median(PETHev, 1, 'omitnan');
end

if summarymeasures(2) == "std"
    PETHd = std(PETHev, 1, 'omitnan');
elseif summarymeasures(2) == "sem"
    PETHd = std(PETHev, 1, 'omitnan')./sqrt(size(PETHev, 1));
elseif summarymeasures(2) == "ci"
    SEM = td(PETHev, 1, 'omitnan')./sqrt(size(PETHev, 1));               % Standard Error
    ts = tinv([0.025  0.975],size(PETHev, 1)-1);      % T-Score
    PETHd = mean(PETHev, 1, 'omitnan') + ts*SEM;                      % Confidence Intervals
elseif summarymeasures(2) == "mad"
    PETHd = mad(PETHev);
end

fig1 = figure(1);
subplot(2, 1, 1)
h = imagesc(PETHts, 1:PETHevnr, (PETHdata));
set(h, 'AlphaData', ~isnan(PETHdata))
% clim([-5 5])
subplot(2, 1, 2)
g = plot(PETHts, median(PETHdata, 1, 'omitnan'), 'LineWidth', 2);
% ylim([-3 3])
xlim([min(PETHts) max(PETHts)])
xline(0, '-r')
% yline(mean(mean(nanBLdata, 1, 'omitnan'), 'omitnan'), '-', 'LineWidth', 1)
titlefig1 = string(strcat(combiname,' Freezing Onset'));
set(fig1, 'Name', titlefig1)
%  title(fig1, string(strcat(combiname,' Freezing Onset')))


fig2 = figure(2);
plot(Time, contvar)
hold on
for ii = 1:size(Events, 1)
    xline(Events(ii, 1), '-g')
    xline(Events(ii, 2), '-r')
end
titlefig2 = string(strcat(combiname,' ZScore with onset events'));
set(fig2, 'Name', titlefig2)



%% 10. Sort freezing episodes by length: 

[~, idx] = sort(boutdur, 'ascend');
%in case of repeated durations:

for  ii= 1:size(boutdur, 1)
    sortedPETH(ii, :) = PETHev((idx(ii)), :);
end

% Criteria for BL (threshold = 0.5 seg):

rowstoelim = underthresh == 1;
if ~isempty(rowstoelim)
    sortedPETH = sortedPETH(~rowstoelim, :);
end


fig3 = figure(3);
tiledlayout(2, 1)
nexttile
h = imagesc(PETHts, 1:PETHevnr, (sortedPETH));
set(h, 'AlphaData', ~isnan(sortedPETH))
% clim([-5 5])
nexttile
plot(PETHts, median(PETHdata, 1, 'omitnan'), 'LineWidth', 2);
xlim([min(PETHts) max(PETHts)])
xline(0, '-r')
% yline(mean(mean(nanBLdata, 1, 'omitnan'), 'omitnan'), '-', 'LineWidth', 1)

% Recortar PETH a rangos deseados: 

idxPETHts = PETHts >= -Pre & PETHts <= Post;
udPETHts = PETHts(idxPETHts);
udPETHdata = PETHdata(:, idxPETHts);

fig4 = figure(4);
tiledlayout(2, 1)
nexttile
h = imagesc(udPETHts, 1:size(udPETHdata, 1), udPETHdata);
set(h, 'AlphaData', ~isnan(udPETHdata))
% clim([-5 5])
nexttile
plot(udPETHts, median(udPETHdata, 1, 'omitnan'), 'LineWidth', 2);
xlim([min(udPETHts) max(udPETHts)])
%ylim([min(udPETHdata) max(udPETHdata)])
xline(0, '-r')
% yline(mean(mean(BLdataZ, 1, 'omitnan'), 'omitnan'), '-', 'LineWidth', 1)


if summarymeasures(1) == "mean"
    udPETHm = mean(udPETHdata, 1, 'omitnan');
elseif summarymeasures(1) == "median"
    udPETHm = median(udPETHdata, 1, 'omitnan');
end

if summarymeasures(2) == "std"
    udPETHd = std(udPETHdata, 1, 'omitnan');
elseif summarymeasures(2) == "sem"
    udPETHd = std(udPETHdata, 1, 'omitnan')./sqrt(size(udPETHdata, 1));
elseif summarymeasures(2) == "ci"
    SEM = td(udPETHdata, 1, 'omitnan')./sqrt(size(udPETHdata, 1));               % Standard Error
    ts = tinv([0.025  0.975],size(udPETHdata, 1)-1);      % T-Score
    udPETHd = mean(udPETHdata, 1, 'omitnan') + ts*SEM;                      % Confidence Intervals
elseif summarymeasures(2) == "mad"
    udPETHd = mad(udPETHdata);
end


%% 11. Mean+std, median+mad,max, min:

% 11.1. Complete data (independent of user):

PETHfsev = bin;

tbl = nan(size(binnedBL, 1), 1); tev = tbl;
meanBL = tbl; stdBL = tbl; medianBL = tbl; madBL = tbl;
maxBL = tbl; minBL = tbl;
meanpersecBL = tbl; stdpersecBL = tbl; medianpersecBL = tbl; madpersecBL = tbl;

meanEV = tbl; stdEV = tbl; medianEV = tbl; madEV = tbl;
maxEV = tbl; minEV = tbl;
meanpersecEV = tbl; stdpersecEV = tbl; medianpersecEV = tbl; madpersecEV = tbl;


for jj = 1:size(binnedBL, 1)
    if double(sum(~isnan(binnedBL(jj, :))))*PETHfsev < 0.5 % criterio de duración.. si ponemos una fsev para el PETH muy alto, hay que tener cuidado porque podría pasar mucho esto al ser un ds
        tbl(jj) = sum(~isnan(binnedBL(jj, :)))*PETHfsev; tev(jj)=sum(~isnan(binnedEV(jj, :)))*PETHfsev;
% BL
        meanBL(jj) = nan; stdBL(jj) = nan; medianBL(jj) = nan; madBL(jj) = nan;
        maxBL(jj) = nan; minBL(jj) = nan;
        meanpersecBL(jj) = nan; stdpersecBL(jj) = nan; medianpersecBL(jj) = nan; madpersecBL(jj) = nan;
% EV
        meanEV(jj) = nan; stdEV(jj) = nan; medianEV(jj) = nan; madEV(jj) = nan;
        maxEV(jj) = nan; minEV(jj) = nan;
        meanpersecEV(jj) = nan; stdpersecEV(jj) = nan; medianpersecEV(jj) = nan; madpersecEV(jj) = nan;        
    else
        tbl(jj) = double(sum(~isnan(binnedBL(jj, :))))*PETHfsev; % duración bin en seg // la dur del bout ya la tenemos
        tev(jj)= double(sum(~isnan(binnedEV(jj, :))))*PETHfsev;
% BL
        meanBL(jj) = mean(binnedBL(jj,:), 'omitnan'); stdBL(jj) = std(binnedBL(jj, :), 'omitnan');
        medianBL(jj) = median(binnedBL(jj,:), 'omitnan'); madBL(jj) = mad(binnedBL(jj, :));
        maxBL(jj) = max(binnedBL(jj,:),[],'omitnan'); minBL(jj) = min(binnedBL(jj,:),[],'omitnan');
        meanpersecBL(jj) = meanBL(jj)/tbl(jj); stdpersecBL(jj) = stdBL(jj)/tbl(jj);
        medianpersecBL(jj) = medianBL(jj)/tbl(jj); madpersecBL(jj) = madBL(jj)/tbl(jj);
% EV
        meanEV(jj) = mean(binnedEV(jj,:), 'omitnan'); stdEV(jj) = std(binnedEV(jj, :), 'omitnan');
        medianEV(jj) = median(binnedEV(jj,:), 'omitnan'); madEV(jj) = mad(binnedEV(jj, :));
        maxEV(jj) = max(binnedEV(jj,:),[],'omitnan'); minEV(jj) = min(binnedEV(jj,:),[],'omitnan');
        meanpersecEV(jj) = meanEV(jj)/tev(jj); stdpersecEV(jj) = stdEV(jj)/tev(jj);
        medianpersecEV(jj) = medianEV(jj)/tev(jj); madpersecEV(jj) = madEV(jj)/tev(jj);
    end    
end


% 11.2. Por intervalos de tiempo determinados por usuario: 

meandata = zeros([length(boutdur) size(AUCint, 1)]); stddata = meandata;
mediandata = meandata; maddata = meandata;maxdata = meandata;mindata = meandata;
durdata = meandata;

for ii = 1:AUCqt
    for jj = 1:size(binnedBL, 1)
        if double(sum(PETHts >= AUCint(ii, 1) & PETHts < AUCint(ii, 2) & ~isnan(PETHev(jj, :))))*bin < 0.5
            meandata(jj, ii) = nan; stddata(jj, ii) = nan;
            mediandata(jj, ii) = nan; maddata(jj, ii) = nan;
            maxdata(jj, ii) = nan; mindata(jj, ii) = nan;
            durdata(jj, ii) = double(sum(~isnan(PETHev(jj, PETHts >= AUCint(ii, 1) & PETHts <= AUCint(ii, 2)))))*bin;
        else
            tempevent = PETHev(jj, (PETHts >= AUCint(ii, 1) & PETHts <= AUCint(ii, 2) & ~isnan(PETHev(jj, :))));
            meandata(jj, ii) = mean(tempevent, 'omitnan'); stddata(jj, ii) = std(tempevent, 'omitnan');
            mediandata(jj, ii) = median(tempevent, 'omitnan'); maddata(jj,ii) = mad(tempevent);
            maxdata(jj, ii) = max(tempevent,[], 'omitnan'); mindata(jj, ii) = min(tempevent,[], 'omitnan');
            durdata(jj, ii) = double(sum(~isnan(PETHev(jj, PETHts >= AUCint(ii, 1) & PETHts <= AUCint(ii, 2)))))*bin;

        end
    end
end  
% 
% sum(PETHts >= AUCint(ii, 1) & PETHts <= AUCint(ii, 2) & ~isnan(PETHev(jj, :)))
% sum(~isnan(PETHev(jj, :)))
% sum(~isnan(PETHev(jj, PETHts >= AUCint(ii, 1) & PETHts <= AUCint(ii, 2))))


% 11. 3. Filtrar los que cumplan la condición que todos tengan la duración
% indicada por usuario. 

intdur = repmat((AUCint(:, 2) - AUCint(:, 1)).', length(boutdur), 1);
indx = double(durdata == intdur); 
indx = sum(indx, 2);

meandatasame = meandata(indx == AUCqt, :); stddatasame = stddata(indx == AUCqt, :); 
mediandatasame = mediandata(indx == AUCqt, :); maddatasame = maddata(indx == AUCqt, :); 
maxdatasame = maxdata(indx == AUCqt, :); mindatasame = mindata(indx == AUCqt, :); 
durdatasame = durdata(indx == AUCqt, :); boutdursame = boutdur(indx == AUCqt, :);


% Tablas para guardar: 

% 1. AUC y dens con su correspondiente duración del completo (BL por un
% lado y Bev por otro)

% 2. AUC y dens con duración de los intervalos seleccionados por usuario
% (independientemente de que el evento sea de menor duración que el
% intervalo establecido)

% 3. AUC y dens de intervalos con filtro de cuales duran lo que tiene que
% durar el intervalo. 

% >>> En todos los datos tengo tantas columnas de AUC, dens y dur como
% intervalos tenga. Además, otra columna con la duración del freezing de
% ese evento. 

% Measurenames: 

% Data 1:

measurenames1 = {'meanpre', 'meanpost','stdpre','stdpost', 'medianpre', 'medianpost',...
 'madpre','madpost', 'maxpre','maxpost', 'minpre', 'minpost','durpre (s)', 'boutdur (s)'};
prepost = table(meanBL, meanEV, stdBL, stdEV, medianBL, medianEV,...
  madBL, madEV, maxBL, maxEV, minBL, minEV, tbl, tev, 'VariableNames', measurenames1 );

% Data 2 and 3 (even if empty):

meannames = {}; stdnames = {}; mediannames = {}; madnames = {};
maxnames = {}; minnames = {}; durnames = {};

for ii = 1:AUCqt % 3 because of three times
    meannames{ii} = sprintf('Mean Int%i', ii); stdnames{ii} = sprintf('Std Int%i', ii);
    mediannames{ii} = sprintf('Median Int%i', ii); madnames{ii} = sprintf('Mad Int%i', ii);
    maxnames{ii} = sprintf('Max Int%i', ii); minnames{ii} = sprintf('Min Int%i', ii);
    durnames{ii} = sprintf('Dur Int%i', ii);
end

% if ~isempty(AUCsamelen)
%     for ii = 1:size(AUCsamelen)
%         AUCsamenames{ii} = sprintf('AUC Int%i', ii);
%         DENSsamenames{ii} = sprintf('DENS Int%i', ii)
%     end
% end
% 

measurenames2 = [meannames stdnames mediannames madnames maxnames minnames durnames 'boutdur'];
allint = array2table([meandata stddata mediandata maddata maxdata mindata durdata boutdur], 'VariableNames', measurenames2);
if ~isempty(meandatasame)
    samelen = array2table([meandatasame stddatasame mediandatasame maddatasame maxdatasame mindatasame durdatasame boutdursame], 'VariableNames', measurenames2);
    %xlsx3 = writetable(samelen);
end

%% Calculate n for each value in mean:

idxPETH = ~isnan(PETHev); nint = sum(idxPETH, 1);
idxudPETH = ~isnan(udPETHdata); nintud = sum(idxudPETH, 1);



%% 12. Save all data and figures: 

% Params:
results.PETH.onset.params.Pre = Pre;
results.PETH.onset.params.Post = Post;
results.PETH.onset.params.Postwfr = Postwfr;
results.PETH.onset.params.PreWind = PreWind;
results.PETH.onset.params.PostWind = PostWind;
results.PETH.onset.params.bin = bin;
results.PETH.onset.params.binsize = binsize; % in samples
results.PETH.onset.params.binthreshold = threshold;
results.PETH.onset.params.summarymeasures = summarymeasures;

% PETH variables: 

results.PETH.onset.PETHdata = PETHev;
results.PETH.onset.PETHwbincriteria = PETHdata;
results.PETH.onset.PETHts = PETHts;
results.PETH.onset.udPETHdata = udPETHdata;
results.PETH.onset.udPETHts = udPETHts;

results.PETH.onset.PETHmean = PETHm;
results.PETH.onset.PETHdev = PETHd;
results.PETH.onset.udPETHmean = udPETHm;
results.PETH.onset.udPETHdev = udPETHd;

results.PETH.onset.sortedPETH = sortedPETH;
results.PETH.onset.EVdata = EVdata;
results.PETH.onset.BLdata = BLdata;

% Figures: 

folderpath1 = strcat(savepath, '/Figures',combiname);
folderpath2 = strcat(savepath, '/Data',combiname);
mkdir(folderpath1, 'onsetPETH');
mkdir(folderpath2, 'onsetPETH')

fig1name = strcat(savepath, '/Figures',combiname,'/onsetPETH/completePETH.jpeg');
fig2name = strcat(savepath, '/Figures', combiname,'/onsetPETH/eventscontvar.jpeg');
fig3name = strcat(savepath, '/Figures',combiname,'/onsetPETH/sortedPETH.jpeg');
fig4name = strcat(savepath, '/Figures',combiname,'/onsetPETH/blthreshPETH.jpeg');

saveas(fig1, fig1name)
saveas(fig2, fig2name)
saveas(fig3, fig3name)
saveas(fig4, fig4name)



% Tables:

xlsx1 = strcat(savepath, '/Data',combiname,'/onsetPETH/prepost.xlsx');
xlsx2 = strcat(savepath, '/Data',combiname,'/onsetPETH/allint.xlsx');
xlsx3 = strcat(savepath, '/Data',combiname,'/onsetPETH/samelengthint.xlsx');
writetable(prepost, xlsx1, 'WriteVariableNames', true); 
writetable(allint, xlsx2, 'WriteVariableNames', true);
if ~isempty(meandatasame)
    writetable(samelen, xlsx3, 'WriteVariableNames', true);
% else
%     message = 'Not enough data from chosen interval (data is not the same length). Data Matrix is empty.';
%     disp(message)
end


% Events:

pethcolcompletets = string(PETHts);
xlsx4 = strcat(savepath, '/Data',combiname,'/onsetPETH/completePETH.xlsx');
xlsx5 = strcat(savepath, '/Data', combiname,'/onsetPETH/sortedPETHwthresh.xlsx');
xlsx6 = strcat(savepath, '/Data',combiname,'/onsetPETH/PETHwthresh.xlsx');
PETHev2 = array2table(PETHev, 'VariableNames', pethcolcompletets);
sortedPETH2 = array2table(sortedPETH, 'VariableNames', pethcolcompletets);
PETHdata2 = array2table(PETHdata, 'VariableNames',pethcolcompletets);

writetable(PETHev2, xlsx4, 'WriteVariableNames', true);
writetable(sortedPETH2, xlsx5, 'WriteVariableNames', true);
writetable(PETHdata2, xlsx6, 'WriteVariableNames', true);

pethcolintts = string(udPETHts);
xlsx7 = strcat(savepath, '/Data',combiname,'/onsetPETH/udPETH.xlsx');
udPETHdata2 = array2table(udPETHdata, 'VariableNames', pethcolintts);
writetable(udPETHdata2, xlsx7, 'WriteVariableNames', true);


% Mean + sem:
summaryPETH = table(PETHts.', PETHm.', PETHd.', nint.', 'VariableNames', {'Time' char(summarymeasures(1)) char(summarymeasures(2)) 'n'});
udsummaryPETH = table(udPETHts.', udPETHm.', udPETHd.', nintud.', 'VariableNames', {'Time' char(summarymeasures(1)) char(summarymeasures(2)) 'n'});

% summaryPETH = table(PETHts.', PETHm.', PETHd.', 'VariableNames', {'Time' char(summarymeasures(1)) char(summarymeasures(2))});
% udsummaryPETH = table(udPETHts.', udPETHm.', udPETHd.', 'VariableNames', {'Time' char(summarymeasures(1)) char(summarymeasures(2))});

xlsx8 = strcat(savepath, '/Data',combiname,'/onsetPETH/summaryPETH.xlsx');
xlsx9 = strcat(savepath, '/Data',combiname,'/onsetPETH/udsummaryPETH.xlsx');

writetable(summaryPETH, xlsx8)
writetable(udsummaryPETH, xlsx9)

close all



end



