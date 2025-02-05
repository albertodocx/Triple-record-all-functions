function bouts = dlcmovement(data,mypath)

if isnan(data)
    disp ('No hay eventos')
    bouts.frames = NaN;
    bouts.seg = NaN;
    bouts.segcut = NaN;

else
    disp('Hay eventos')

    separacion_minima = 15;
    mypath = strcat(mypath,'/conducta');

    % Calcular la envolvente del movimiento (suavizado)
    smoothed_data = smooth(abs(data),100);
    smoothed_data = rescale(smoothed_data,0,1);
    datares = rescale(data,0,1);
    umbral = 0.09;

    % Encontrar los momentos donde la magnitud de la derivada supera el umbral
    movimientos = abs(smoothed_data) > umbral;

    % Determinar los índices donde se produce el movimiento
    indices_movimiento = find(movimientos);

    % Inicializar los intervalos de movimiento
    intervalos_movimiento = [];

    % Identificar los intervalos continuos de movimiento
    inicio_intervalo = indices_movimiento(1);
    for i = 2:length(indices_movimiento)
        if indices_movimiento(i) - indices_movimiento(i-1) > separacion_minima
            % Se alcanzó el final de un intervalo de movimiento
            fin_intervalo = indices_movimiento(i-1);
            intervalos_movimiento = [intervalos_movimiento; inicio_intervalo, fin_intervalo];
            % Iniciar un nuevo intervalo de movimiento
            inicio_intervalo = indices_movimiento(i);
        end
    end
    % Agregar el último intervalo de movimiento
    fin_intervalo = indices_movimiento(end);
    intervalos_movimiento = [intervalos_movimiento; inicio_intervalo, fin_intervalo];
    
    min_ev = min(data);
    max_ev = max(data);
    max_range = max(abs(min_ev), abs(max_ev));
    y_limit = [0, max_range];
    graph_limit = 1;%max(smoothed_data);

    % Hacer la figura
    figure('Position',[111,334,1227,420]);

    subplot(2,1,1)
    plot(datares,"Color",[0.97,0.82,1.00]);
    hold on;
    plot(smoothed_data,"Color",[0.35,0.13,0.39],"LineWidth",1.5);
    hold on;
    xdddlim = [1, numel(data)];
    line(xdddlim, [umbral, umbral], 'Color', 'k', 'LineStyle', '-',"LineWidth",1);
    hold off

    title('Detección de intervalos de movimiento');
    ylim([0, graph_limit])
    xlabel('Tiempo');
    ylabel('Valor de movimiento');
    legend('Datos brutos','Envolvente','Umbral');

    subplot(2,1,2)
    for i = 1:size(intervalos_movimiento, 1)
        rectangle('Position', [intervalos_movimiento(i, 1), y_limit(1), intervalos_movimiento(i, 2) - intervalos_movimiento(i, 1),diff(y_limit)], 'FaceColor', [0.8,0.8,0.8], 'EdgeColor', 'none');
        hold on;
    end

    plot(smoothed_data,"Color",[0.35,0.13,0.39],"LineWidth",1.5);
    hold off

    title('Intervalos de movimiento');
    ylim([0, graph_limit])
    xlabel('Tiempo');
    ylabel('Valor del movimiento');
    legend('Envolvente');
    nombre_grafico = fullfile(mypath, 'boutsconducta.png');
    saveas(gcf, nombre_grafico);


    segundos = intervalos_movimiento/15;
    resta = segundos(:, 2) - segundos(:, 1);
    filas_validas = resta >= 1;
    segundoscut = segundos(filas_validas, :);

    % Devolver los intervalos de movimiento como resultado
    bouts.frames = intervalos_movimiento;
    bouts.seg = segundos;
    bouts.segcut = segundoscut;
    bouts.data.raw = data;
    bouts.data.smoothed = smoothed_data;
    bouts.data.rescaled = datares;
end
end
%
%     % Encontrar los picos de la envolvente
%     [~, locs] = findpeaks(smoothed_data, 'MinPeakHeight', umbral, 'MinPeakDistance', separacion_minima);
%
%     % Identificar los intervalos continuos de movimiento
%     intervalos_movimiento = [];
%     inicio_intervalo = locs(1);
%
%     for i = 2:length(locs)
%         if locs(i) - locs(i-1) > separacion_minima
%             % Se alcanzó el final de un intervalo de movimiento
%             fin_intervalo = locs(i-1);
%             intervalos_movimiento = [intervalos_movimiento; inicio_intervalo, fin_intervalo];
%             % Iniciar un nuevo intervalo de movimiento
%             inicio_intervalo = locs(i);
%         end
%     end
%
%     % Agregar el último intervalo de movimiento
%     fin_intervalo = locs(end);
%     intervalos_movimiento = [intervalos_movimiento; inicio_intervalo, fin_intervalo];
%
%     % Graficar la señal suavizada y los intervalos de movimiento
%     figure;
%     plot(abs(data));
%     hold on;
%     plot(smoothed_data, 'r', 'LineWidth', 1.5);
%     hold on;
%     for i = 1:size(intervalos_movimiento, 1)
%         plot(intervalos_movimiento(i,:), smoothed_data(intervalos_movimiento(i,:)), 'g', 'LineWidth', 2);
%     end
%     title('Detección de intervalos de movimiento con envolvente');
%     xlabel('Índice');
%     ylabel('Valor de la envolvente');
%     legend('Señal original', 'Envolvente', 'Intervalos de Movimiento');
%     nombre_grafico = fullfile(mypath, 'boutsconducta.png');
%     saveas(gcf, nombre_grafico);
%
%     % Convertir los índices a segundos
%     segundos = locs / 15;
%     resta = segundos(2:end) - segundos(1:end-1);
%
%     % Filtrar intervalos que tengan una duración mayor o igual a 1 segundo
%     filas_validas = resta >= 1;
%     segundoscut = [segundos(filas_validas), segundos([false; filas_validas])];
%
%     % Devolver los intervalos de movimiento como resultado
%     bouts.frames = locs;
%     bouts.seg = segundos;
%     bouts.segcut = segundoscut;
% end
% end

