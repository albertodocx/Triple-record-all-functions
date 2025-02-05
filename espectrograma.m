function espectrograma(TDTdata,results,mypath)

% imagesc(1:1729,psprofilebueno.pSfreqs,pow2db(psprofilebueno.pS))
% ylim([100 1000]);
% clim([0 70]);
% colormap("jet")

% Obtener información de los datos
fs = results.lfp.psprofilebueno.pSfreqs;
T=[1 results.lfp.psprofilebueno.params.nWin];
PS=pow2db(results.lfp.psprofilebueno.pS);
    % % Calculate spectrogram
    % [S, F, T] = spectrogram(datosFiltrados, hann(2), 1, 1, fs);
    % 
    % % Convert to decibels using pow2db
    % S_dB = max(0, pow2db(abs(S))); % Limit to 0 dB or a minimum value if desired

    % Plot spectrogram
    figure;
    imagesc(T, results.lfp.psprofilebueno.pSfreqs, PS);
    title('Spectrogram (dB)');
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    colorbar;
    colormap('jet');
    set(gca, 'YDir', 'normal');
    pbaspect([1 0.5 1]);
    clim([0 100]);

    % Save the spectrogram
    if nargin > 2
        saveas(gcf, fullfile(mypath, 'spectrogram.png'));
    end

    % Create and plot another spectrogram with a frequency range of 0.5 to 15 Hz
    figure;
    freqRange = (fs >= 0.5) & (fs <= 15);
    imagesc(T, fs(freqRange), PS(freqRange, :));
    title('Spectrogram (0.5-15 Hz) (dB)');
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    colorbar;
    colormap('jet');
    set(gca, 'YDir', 'normal');
    pbaspect([1 0.5 1]);
    clim([0 100]);

    % Save the second spectrogram
    if nargin > 2
        saveas(gcf, fullfile(mypath, 'spectrogram_0.5-15Hz.png'));
    end
end


% function espectrograma(TDTdata, datosFiltrados, mypath)
%     % Obtener información de los datos
%     fs = TDTdata.streams.Wav1.fs;
% 
%     % Calculate spectrogram
%     [S, F, T] = spectrogram(datosFiltrados, hann(256), 128, 1024, fs);
% 
%     % Plot spectrogram
%     figure;
%     imagesc(T, F, 10*log10(abs(S)));
%     title('Spectrogram');
%     xlabel('Time (s)');
%     ylabel('Frequency (Hz)');
%     colorbar;
%     colormap('jet');
%     set(gca, 'YDir', 'normal');
%     pbaspect([1 0.5 1]);
% 
%     % Save the spectrogram
%     if nargin > 2
%         saveas(gcf, fullfile(mypath, 'spectrogram.png'));
%     end
% 
%     % Create and plot another spectrogram with a frequency range of 0.5 to 15 Hz
%     figure;
%     freqRange = (F >= 0.5) & (F <= 15);
%     imagesc(T, F(freqRange), 10*log10(abs(S(freqRange, :))));
%     title('Spectrogram (0.5-15 Hz)');
%     xlabel('Time (s)');
%     ylabel('Frequency (Hz)');
%     colorbar;
%     colormap('jet');
%     set(gca, 'YDir', 'normal');
%     pbaspect([1 0.5 1]);
% 
%     % Save the second spectrogram
%     if nargin > 2
%         saveas(gcf, fullfile(mypath, 'spectrogram_0.5-15Hz.png'));
%     end
% end
