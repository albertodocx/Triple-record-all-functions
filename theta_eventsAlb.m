function thetaanalisis = theta_eventsAlb(results,mypath)

    thetach_profile = results.lfp.psprofilebueno.pS;
    freqs = results.lfp.psprofilebueno.pSfreqs;
    
    %theta_cf = results.psProfile.params.theta_cf;
    theta_cf = [5.5 8.5]; % Mismo rango que Alberto
    thetafq = (freqs >= theta_cf(1) & freqs <= theta_cf(2)); % fq range
    thetarg = freqs(thetafq);
    thetach_psth = thetach_profile(thetafq, :); % Espectro sólo para las fq de interés
    
    %delta_cf = results.psProfile.params.delta_cf;
    delta_cf = [1 3.5]; % No puede ser de menos por la ventana de desplazamiento que hemos puesto en power spectrum (2s)
    deltafq = (freqs >= delta_cf(1) & freqs <= delta_cf(2)); % based on Malezieux, M., Kees, A. L., & Mulle, C. (2020)
    deltarg = freqs(deltafq);
    thetach_psde= thetach_profile(deltafq, :);

    if size(thetach_psde, 1) > size(thetach_psth, 1)
        thetach_psde = thetach_psde(1:size(thetach_psth, 1), :);
    elseif size(thetach_psde, 1) < size(thetach_psth, 1)
        thetach_psth = thetach_psth(1:size(thetach_psde, 1), :);
    end
    
    time_spec = (0:size(thetach_psth, 2)-1)./2; % 2 muestras por seg

% Ratio theta/delta
    thetamean = mean(thetach_psth, 1); % Vector de tiempo >> media de fq
    deltamean = mean(thetach_psde, 1);
    ratiotd = thetamean./deltamean; 
    
    % ratiotd = deltamean./thetamean;
    if sum(isinf(ratiotd)) > 0
        ratiotd(ratiotd == Inf) = nan;
    end

    % Código de Alberto para sgolayfilt:
    ratiofilt = movmean(movmean( sgolayfilt(ratiotd,50,51), 20), 20); % Sale igual que si lo normalizas antes
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

%Detect run-theta events
    LIAthresh = 0.2; % save to detect ripples events under this threshold
    runthresh = 0.3; % Estos valores después de normalizar entre 0 y 1

    % Compare SP and SLM to check which channel to choose for theta calc.: 

    imrg = freqs >= 1 & freqs <= 15;
    imgfreqs = freqs(imrg);
    
    fig2 = figure (2);
    tiledlayout(3, 1)
    
    %plot spectrogram
    nexttile
    imagesc(time_spec, imgfreqs, pow2db(abs(thetach_profile(imrg, :))))
    title('Power Spectrum (dB)')
    ylabel('Fq (Hz)')
    xlabel('Time (s)')
    colorbar
    colormap('jet')
    set(gca, 'YDir', 'normal')
    clim([0 100])

    %plot normalized theta power
    nexttile % Plot normalized theta power
    plot(time_spec, thetameansc, 'Color' ,'#dad0d8')
    hold on
    plot(time_spec, thetafilt, 'blue')
    title('Normalized Theta Power')
    xlabel('Time (s)')
    hold off

    %plot normalized ratio power
    nexttile % Plot normalized ratio power
    plot(time_spec, ratiosc, 'Color' ,'#dad0d8')
    hold on
    plot(time_spec, ratiofilt, 'blue')
    title('Normalized Ratio Theta/Delta Power')
    xlabel('Time (s)')
    yline(LIAthresh, 'red')
    yline(runthresh, 'green')
    hold off
    
    fig2name = strcat(mypath, '/specthtratio.png');
    saveas(fig2, fig2name);
    close all

    fig3 = figure (3);
    plot(time_spec, thetameansc, 'Color' ,'#dad0d8')
    hold on
    plot(time_spec, thetafilt, 'blue')
    xlabel('Time (s)')
    xlim ([0, 300])
    ylim([0, 1])
    pbaspect([1 0.5 1]);
    hold off

    fig3name = strcat(mypath, '/normtheta.png');
    saveas(fig3, fig3name);
    close all

    fig4 = figure (4);
    plot(time_spec, ratiosc, 'Color' ,'#dad0d8')
    hold on
    plot(time_spec, ratiofilt, 'blue')
    xlabel('Time (s)')
    yline(LIAthresh, 'red')
    yline(runthresh, 'green')
    xlim ([0, 300])
    ylim([0, 1])
    pbaspect([1 0.5 1]);
    hold off

    fig4name = strcat(mypath, '/ratio.png');
    saveas(fig4, fig4name);
    close all

    fig5 = figure (5);
    nexttile
    imagesc(time_spec, imgfreqs, pow2db(abs(thetach_profile(imrg, :))))
    ylabel('Fq (Hz)')
    xlabel('Time (s)')
    colorbar
    colormap('jet')
    set(gca, 'YDir', 'normal')
    clim([0 100])
    xlim ([0, 300])
    pbaspect([1 0.5 1]);

    fig5name = strcat(mypath, '/spectrum.png');
    saveas(fig5, fig5name);
    close all

% Get number of events above that threshold:
 
    thetabouts = ratiofilt > runthresh;
        if thetabouts(1) == 1
            [~, locs] = findpeaks(double(thetabouts)); % Para calcular el núm de eventos de theta (realmente te coge el inicio de cada uno de ellos).
            locs = [1 locs];
        else
            [~, locs] = findpeaks(double(thetabouts)); % Para calcular el núm de eventos de theta (realmente te coge el inicio de cada uno de ellos).
        end
       
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

    boutsseg = bouts/2;
    
    % Determina el tiempo máximo
    tiempo_maximo = max(time_spec);

% Crea un gráfico de eventos usando rectangles
    figure ('Name','Event: on/off');
    
    for i = 1:size(boutsseg, 1)
        rectangle('Position', [boutsseg(i, 1), 0, boutsseg(i, 2) - boutsseg(i, 1), 1], 'FaceColor', 'k', 'EdgeColor', 'none');
        hold on;
    end
    
    title('Gráfico de Eventos');
    xlabel('Tiempo');
    ylabel('Estado del evento');
    xlim([0, tiempo_maximo]); % Establece el límite del tiempo máximo
    ylim([0, 1]); % Establece límites para un mejor aspecto visual
    grid on;
    saveas(gcf, fullfile(mypath, 'events.png'));
    
    figure ('Name','Event: on/off');
    
    for i = 1:size(boutsseg, 1)
        rectangle('Position', [boutsseg(i, 1), 0, boutsseg(i, 2) - boutsseg(i, 1), 1], 'FaceColor', 'k', 'EdgeColor', 'none');
        hold on;
    end
    
    xlabel('Tiempo');
    ylabel('Estado del evento');
    xlim([0, 300]); % Establece el límite del tiempo máximo
    ylim([0, 1]); % Establece límites para un mejor aspecto visual
    pbaspect([1 0.5 1]);
    grid on;
    saveas(gcf, fullfile(mypath, 'eventscut.png'));
    

%Guardar en un struct
    thetaanalisis.profile = thetach_profile;
    thetaanalisis.freqs = freqs;
    thetaanalisis.bouts.muestras = bouts;
    thetaanalisis.bouts.segundos = boutsseg;
    thetaanalisis.ratio.params.LIA = LIAthresh;
    thetaanalisis.ratio.params.RUN = runthresh;
    thetaanalisis.ratio.thetadelta = ratiosc;

    
end