

ripplech = results.params.ripplech;
LIAthresh = results.thetaevent.params.LIAthresh;
rippleps = results.psProfile.pS(:, :, ripplech);

pSfreqs = results.psProfile.pSfreqs;

spratio = results.thetaevent.spthetafilt; 

plot(spratio)
hold on
yline(LIAthresh, 'red')
hold off

LIAevents = spratio <= LIAthresh;


% Filtrar datos de señal downsampleada de canal de ripples: 

LFPripplech = results.LFPds.dsdata(:, ripplech);
Fs = results.LFPds.fs;
LFPds = results.LFPds.dsdata;

pass1 = 90; % Hz
pass2 = 250; % Hz
%order = 1700; % como en la función filtroFIR1 (por si diera error)

[signRip_complete, signEnv_complete] = filtroFIR1(Fs, LFPds, pass1, pass2);
[signRip_ripples, signEnv_ripples] = filtroFIR1(Fs, LFPripplech, pass1, pass2);

plot(signRip_ripples(1000000:1008000))
hold on 
plot(signEnv_ripples(1000000:1008000))
hold off


% Generate empty vector

filtripplesmean = mean(signRip_ripples);

wLIAthresh = repelem(filtripplesmean, size(LFPds, 1));

wLIAthresh(LIAevents) = signRip_ripples(LIAevents);

plot(wLIAthresh(1:8000))

LIAevents(1000000:1008000)

