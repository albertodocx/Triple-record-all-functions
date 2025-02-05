function [ pS, band, freq, h, h2 ] = powerSpectrumProfileAlb( d, Fs, timeSpec )
% Power channel profile and power spectrum over 2s windows + hanning
% [pS, band] = powerSpectrumProfile( d, Fs )
% INPUT
%   d: data matrix with all sorted channel 
%   Fs: sampling frequency (scalar)
%   timeSpec: do time spectrum in max theta channel 0 (default) or 1 (do it)
% OUTPUT
%   pS: 3d power spectrum matrix (frequencies x time windows (2s) x channel)
%   band: structure with frequencies band per channel(delta, theta, gamma, hfos).
%   freq: frequencies vector.
%   h: profile image handle 
%   h2: power spectrum image handle
% LCN-MV 2016
% Edit: LCN-AN 2019

% Electrode position, mm with respect bregma

% Data matrix needs to be time x channels
if size(d,2) > size(d,1) 
    d = d';
end
% Number of channels
nCh = size(d,2);

% Default: timeSpec = 0
if ~exist('timeSpec') | isempty(timeSpec)
    timeSpec = 0;
end

% Normalization to zero mean for each channel
d = d - mean(d,1);


% DO POWER SPECTRUM IN WINDOWS OF 2s EVERY 0.5s FOR EACH CHANNEL
% 2 second window
tWin = Fs*2;
% Window moving time step
dtWin = tWin / 4;
% Number of windows
nWin = floor(size(d,1)/dtWin) - 3;
% Power Spectrum storing matrix

 pS = zeros( Fs-1, nWin, nCh );
%pS = struct();

% For each channel
for ch = 1:nCh % Recorro todos los canales.
    for win = 1:nWin
        % Temporal data of the corresponding window
        dataWin = d( dtWin*(win-1)+1 : dtWin*(win-1)+tWin , ch );
        % Storing results of power spectrum between 1 and 1000
        [ pS( :, win, ch), freq ] = power_spectrum( dataWin', true, Fs, [1 1000]);
%         [ pS(:, win, ch), freq ] = pspectrum( dataWin.',Fs,'Reassign', true, 'FrequencyLimits', [1 1000]);

    end
end

% GET POWER OF MAIN RHYTHMS: THETA, GAMMA, HFOs FOR EACH WINDOW AND EACH CHANNEL
% For each channel
for ch = 1:nCh 
    % For each window
    for win = 1:nWin
        % Delta power: between 0 and 3Hz
         pdelta(win,ch) = pow2db( sum( pS( freq>=0 & freq<=3, win, ch )));

        % Theta power: between 4 and 12Hz
         ptheta(win,ch) = pow2db( sum( pS( freq>=4 & freq<=12, win, ch )));

        % Alpha power: between 9 and 16Hz
         palpha(win,ch) = pow2db( sum( pS( freq>=9 & freq<=16, win, ch )));

        % Gamma power: between 20 and 90Hz
        pgamma(win,ch)  = pow2db( sum( pS( freq>=20 & freq<=90, win, ch )));

        % HFO power: between 100 and 1000Hz
        try phfo(win,ch) = 10*log10( sum( pS( freq>=90 & freq<=250, win, ch ))); catch phfo = []; end

    end
end
% Means and stds for each rhythm
meanDelta = mean(pdelta);       stdDelta = std(pdelta);     chDelta = find(meanDelta ~= -Inf );
meanTheta = mean(ptheta);       stdTheta = std(ptheta);     chTheta = find(meanTheta ~= -Inf);
meanAlpha = mean(palpha);       stdAlpha = std(palpha);     chAlpha = find(meanAlpha ~= -Inf );
meanGamma = mean(pgamma);       stdGamma = std(pgamma);     chGamma = find(meanGamma ~= -Inf);
meanHfo   = mean(phfo);         stdHfo   = std(phfo);       chHfo   = find(meanHfo ~= -Inf);


% MAKE FIGURE
h = figure; hold on;
% Plot std with shadows
i1a=fill([chHfo flip(chHfo) chHfo(1)],...
    [meanHfo(chHfo)+stdHfo(chHfo) flip(meanHfo(chHfo)-stdHfo(chHfo))...
    meanHfo(chHfo(1))+stdHfo(chHfo(1))],[.1 .1 .1],'facealpha',.1,'lineStyle','none');
i2a=fill([chGamma flip(chGamma) chGamma(1)],...
    [meanGamma(chGamma)+stdGamma(chGamma) flip(meanGamma(chGamma)-stdGamma(chGamma))...
    meanGamma(chGamma(1))+stdGamma(chGamma(1))],[.2 .9 .2],'facealpha',.1,'lineStyle','none');
i3a=fill([chTheta flip(chTheta) (chTheta(1))],...
    [meanTheta(chTheta)+stdTheta(chTheta) flip(meanTheta(chTheta)-stdTheta(chTheta))...
    meanTheta(chTheta(1))+stdTheta(chTheta(1))],[.7 .2 .2],'facealpha',.1,'lineStyle','none');
i4a=fill([chDelta flip(chDelta) (chDelta(1))],...
    [meanDelta(chDelta)+stdDelta(chDelta) flip(meanDelta(chDelta)-stdDelta(chDelta))...
    meanDelta(chDelta(1))+stdDelta(chDelta(1))],[.2 .2 .7],'facealpha',.1,'lineStyle','none');
i5a=fill([chAlpha flip(chAlpha) (chAlpha(1))],...
    [meanAlpha(chAlpha)+stdAlpha(chAlpha) flip(meanAlpha(chAlpha)-stdAlpha(chAlpha))...
    meanAlpha(chAlpha(1))+stdAlpha(chAlpha(1))],[.7 .2 .7],'facealpha',.1,'lineStyle','none');
% Plot means
i1 = plot( chHfo, meanHfo(chHfo), '-o', 'color', [.1 .1 .1], 'LineWidth' , 2, 'Markersize', 8, 'MarkerFaceColor', 'w' );
i2 = plot( chGamma, meanGamma(chGamma), '-o', 'color', [.2 .9 .2], 'LineWidth' , 2, 'Markersize', 8, 'MarkerFaceColor', 'w' );
i3 = plot( chTheta, meanTheta(chTheta), '-o', 'color', [.7 .2 .2], 'LineWidth' , 2, 'Markersize', 8, 'MarkerFaceColor', 'w' );
i4 = plot( chDelta, meanDelta(chDelta), '-o', 'color', [.2 .2 .7], 'LineWidth' , 2, 'Markersize', 8, 'MarkerFaceColor', 'w' );
i5 = plot( chAlpha, meanAlpha(chAlpha), '-o', 'color', [.7 .2 .7], 'LineWidth' , 2, 'Markersize', 8, 'MarkerFaceColor', 'w' );
% Legend
ylabel('Power (dB)'); xlabel('# Channel');
ax = axis; xlim([1 ax(2)]);
legend([i1 i2 i3 i4 i5],'HFOs 100-500 Hz','Gamma 20-90 Hz','Theta 4-12 Hz',...
    'Alpha 9-15 Hz','Delta 1-3 Hz','location','southeast');


% IF timeSpec IS TRUE, PLOTS THE POWER SPECTRUM OF THETA
if timeSpec
    % Computes channel with maximum Theta
    fprintf('\n -> Channel with max theta is %d .\n',(find(max(meanTheta)==meanTheta)));
    % Question: do you want to see the power Spectrum alon time for the channel that has more theta power?
    qu = input('Power spectrum of this channel (y/n)? ','s');
    switch qu
        % If you think is correct, then 'yes', to distinguish between theta and LIA
        case 'y'
            chSLM = find(max(meanTheta)==meanTheta);
        % If you think SLM is another channel, then 'no'
        case 'n'
            chSLM = input('channel for power spectrum?: ');
        otherwise
            chSLM = find(max(meanTheta)==meanTheta);
    end
    hold off

    % En cada tags calculo el power de theta (columna 3 de tags) y de delta (columna 4 de tags)
    % y hago el cociente (Senior et al., 2008) (columna 5 de tags)
    % Power spectrum of SLM in decibels
%     pSslm = 10 * log10( pS(:,:,chSLM) ) + 60;
    pSslm = pow2db( pS(:,:,chSLM) );

    % Compute relation between theta and LIA powers (Senior et al., 2008)
    tvsd = [];
    for win = 1:nWin
        tvsd(win,1) = sum( pSslm( find( freq>=5 & freq<=9), win ) );
        tvsd(win,2) = sum( pSslm( find( freq>=1 & freq<=3), win ) );
        tvsd(win,3) = tvsd(win,1) / tvsd(win,2);
    end

    h2 = figure; % Hago un espectograma y le ploteo encima la relaci�n theta/delta
    subplot(7,1,2:7)
    h1 = imagesc([1 length(d)/Fs],(freq(find(freq>=1 & freq<=30))),(pSslm((find(freq>=1 & freq<=30)),:)));
    xlabel('s)');
    ylabel('freq(Hz)');
    colormap jet
    h11 = colorbar('southoutside'); h11.Label.String = 'Power (dB)';
    subplot(7,1,1)
    plot(tvsd(:,3));
    xlim([1 length(tvsd)]); ylim([1 3]);
    ylabel('th/del');
end

% Saving info
% Theta
band.theta.pmean = meanTheta;
band.theta.pstd = stdTheta;
band.theta.ch = chTheta;
% Delta
band.delta.pmean = meanDelta;
band.delta.pstd = stdDelta;
band.delta.ch = chDelta;
% Gamma
band.gamma.pmean = meanGamma;
band.gamma.pstd = stdGamma;
band.gamma.ch = chGamma;
% HFO
band.hfo.pmean = meanHfo;
band.hfo.pstd = stdHfo;
band.hfo.ch = chHfo;
% Alpha
band.alpha.pmean = meanAlpha;
band.alpha.pstd = stdAlpha;
band.alpha.ch = chAlpha;

end


