function datosFiltrados = notchfilter(TDTdata)
    % Acceder a los datos
    electrofisiologia = TDTdata.currentsignals.Electro;
    fs = TDTdata.streams.Wav1.fs;    

    % Especificar los parámetros del filtro notch
    wo = 50 / (fs/2);
    factor = 10;
    bw = wo / factor;
    
    % Diseñar el filtro notch
    [b, a] = iirnotch(wo, bw);
    
    % Aplicar el filtro notch a los datos
    datosFiltrados = filtfilt(b, a, electrofisiologia);
    
    % Graficar los datos originales y los datos filtrados
    figure;
    plot(electrofisiologia, 'b', 'LineWidth', 1.5);
    hold on;
    plot(datosFiltrados, 'r', 'LineWidth', 1.5);
    title('Datos Originales y Datos Filtrados');
    legend('Datos Originales', 'Datos Filtrados');
    xlabel('Muestras');
    ylabel('Amplitud');
    grid on;
    
end

