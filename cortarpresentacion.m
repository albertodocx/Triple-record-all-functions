function cortarpresentacion(results,mypath)
    
    % Acceder a los datos de tiempo y fotometría
    tiempo = results.FP.Signals.raw.Time;
    fotometria = results.FP.Signals.DFFModZscore;

    % Valor objetivo
    valor_objetivo = 300;

    % Especificar la tolerancia
    tolerancia = 1e-3;

    % Obtener el número de columnas
    columnas = length(tiempo);

    % Inicializar variables
    columna_encontrada = [];
    distancia_minima = inf;

    % Iterar sobre las columnas para encontrar la más cercana
    for col = 1:columnas
        distancia_actual = abs(tiempo(col) - valor_objetivo);

        if distancia_actual < distancia_minima && distancia_actual < tolerancia
            columna_encontrada = col;
            distancia_minima = distancia_actual;
        end
    end

    % Mostrar el resultado
    if ~isempty(columna_encontrada)
        disp(['El valor más cercano a ', num2str(valor_objetivo), ' se encuentra en la columna ', num2str(columna_encontrada)]);

        % Quedarse con los datos hasta la columna encontrada
        tiempo_cortado = tiempo(1:columna_encontrada);
        fotometria_cortada = fotometria(1:columna_encontrada);

        % Crear una tabla con los datos
        tabla_datos = table(tiempo_cortado', fotometria_cortada', 'VariableNames', {'Tiempo', 'Fotometria'});

        % Replicar el valor de results.FP.path para que coincida con la longitud de la tabla
        valor_path = results.FP.path;
        tabla_datos.Path = repmat({valor_path}, size(tabla_datos, 1), 1);

        % Obtener el nombre del archivo Excel
        nombre_archivo_excel = fullfile(mypath, 'resultadosfotometríacut.xlsx');

        % Guardar la tabla en un archivo Excel
        writetable(tabla_datos, nombre_archivo_excel);

        disp(['Los resultados se han guardado en: ', nombre_archivo_excel]);

        % Crear un gráfico
        figure;
        
        % Ajustar el rango en el eje x
        rango_x = [0 300];
        rango_y = [-4 5];
        plot(tiempo_cortado, fotometria_cortada);
        xlim(rango_x);
        ylim (rango_y);
        pbaspect([1 0.5 1])
        
        xlabel('Tiempo');
        ylabel('Fotometría');
               
        % Guardar el gráfico en el mismo archivo Excel
        nombre_grafico = fullfile(mypath, 'fotometríacut.png');
        saveas(gcf, nombre_grafico);
        disp(['El gráfico se ha guardado en: ', nombre_grafico]);

    else
        disp(['No se encontró ningún valor cercano a ', num2str(valor_objetivo)]);
        
        % Crear una tabla con los datos
        tabla_datos = table(tiempo', fotometria', 'VariableNames', {'Tiempo', 'Fotometria'});

        % Replicar el valor de results.FP.path para que coincida con la longitud de la tabla
        valor_path = results.FP.path;
        tabla_datos.Path = repmat({valor_path}, size(tabla_datos, 1), 1);

        % Obtener el nombre del archivo Excel
        nombre_archivo_excel = fullfile(mypath, 'resultadosfotometríacut.xlsx');

        % Guardar la tabla en un archivo Excel
        writetable(tabla_datos, nombre_archivo_excel);

        disp(['Los resultados se han guardado en: ', nombre_archivo_excel]);

        % Crear un gráfico
        figure;
        
        % Ajustar el rango en el eje x
        rango_x = [0 300];
        rango_y = [-4 5];
        plot(tiempo, fotometria);
        xlim(rango_x);
        ylim (rango_y);
        pbaspect([1 0.5 1])

        xlabel('Tiempo');
        ylabel('Fotometría');
        
        % Guardar el gráfico en el mismo archivo Excel
        nombre_grafico = fullfile(mypath, 'fotometríacut.png');
        saveas(gcf, nombre_grafico);
        disp(['El gráfico se ha guardado en: ', nombre_grafico]);
            end
end


