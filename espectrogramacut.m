function espectrogramacut(TDTdata, datosFiltrados, mypath)
    % Obtener información de los datos
    fs = TDTdata.streams.Wav1.fs;

    % Duración deseada en segundos
    desiredDuration = 300;

    % Número de puntos que corresponden a la duración deseada
    desiredPoints = round(desiredDuration * fs);

    % Asegurarse de que no exceda la longitud de los datos
    if desiredPoints > length(datosFiltrados)
        % Rellenar con ceros si los datos son más cortos que la duración deseada
        datosFiltrados(end+1:desiredPoints) = 0;
    else
        % Tomar solo los primeros desiredPoints de datos
        datosFiltrados = datosFiltrados(1:desiredPoints);
    end

    % Calcular el espectrograma
    [S, F, T] = spectrogram(datosFiltrados, hann(256), 128, 1024, fs);

    % Plot spectrogram
    figure;
    imagesc(T, F, 10*log10(abs(S)));
   
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    colorbar;
    colormap('jet');
    set(gca, 'YDir', 'normal');
    pbaspect([1 0.5 1]);

    % Establecer manualmente los límites de la escala de colores
    clim([-17 22]);

    % Save the spectrogram
    if nargin > 2
        saveas(gcf, fullfile(mypath, 'spectrogram_cut.png'));
    end

    % Crear y representar otro espectrograma con un rango de frecuencia de 0.5 a 15 Hz
    figure;
    freqRange = (F >= 0.5) & (F <= 15);
    imagesc(T, F(freqRange), 10*log10(abs(S(freqRange, :))));
    
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    colorbar;
    colormap('jet');
    set(gca, 'YDir', 'normal');
    pbaspect([1 0.5 1]);

    % Establecer manualmente los límites de la escala de colores
    clim([-17 22]);

    % Save the second spectrogram
    if nargin > 2
        saveas(gcf, fullfile(mypath, 'spectrogram_0.5-15Hz_cut.png'));
    end
end


