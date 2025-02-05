function DLC = openDLC(mypath)
    % Interactuar con el usuario para ingresar nombres de columnas
    aaaprompt = {'bodypart 1', 'bodypart 2', 'bodypart 3', 'object A'};
    % dlg_title = 'Ingrese los nombres base para las categorías';
    % num_lines = 1;
    defaultans = {'Mano derecha', 'Mano izquierda', 'Cola','Pie'};
    base_names = defaultans;%inputdlg(prompt, dlg_title, num_lines, defaultans);
    names = vertcat(aaaprompt,base_names);
    

    % Generar nombres de columnas basados en las bases ingresadas
    column_names = cell(1, numel(base_names) * 3);
    for i = 1:numel(base_names)
        column_names{3*(i-1) + 1} = ['X ' base_names{i}];
        column_names{3*(i-1) + 2} = ['Y ' base_names{i}];
        column_names{3*(i-1) + 3} = ['prob ' base_names{i}];
    end 

    titulos = ['Frames',column_names];

    % Obtiene la lista de archivos en la carpeta
    archivos = dir(fullfile(mypath, '*.csv'));

    % Crea una carpeta para guardar las gráficas si no existe
    if ~exist(fullfile(mypath, 'conducta'), 'dir')
        mkdir(fullfile(mypath, 'conducta'));
    end

    % Verifica si se encontró exactamente un archivo CSV
    if numel(archivos) == 1
        % Lee el archivo CSV como una tabla y asigna los títulos de las columnas
        datos = readtable(fullfile(mypath, archivos(1).name), 'ReadVariableNames', false);
        datos.Properties.VariableNames = titulos;

        % Extrae las coordenadas de los puntos para cada objeto
        x_manodch = datos{:, 2};
        y_manodch = datos{:, 3};
        x_manoizq = datos{:, 5};
        y_manoizq = datos{:, 6};
        x_cola = datos{:, 8};
        y_cola = datos{:, 9};
        x_pie = datos{:, 11};
        y_pie = datos{:, 12};

        % Extrae las probabilidades de los puntos para cada objeto
        prob_manodch = datos{:, 4};
        prob_manoizq = datos{:, 7};
        prob_cola = datos{:, 10};
        prob_pie = datos{:, 13};

        % Grafica los puntos
        figure;
        hold on;

        scatter(x_manodch(prob_manodch>0.9), -y_manodch(prob_manodch>0.9), 'b', 'filled');
        scatter(x_manoizq(prob_manoizq>0.9), -y_manoizq(prob_manoizq>0.9), 'r', 'filled');
        scatter(x_cola(prob_cola>0.9), -y_cola(prob_cola>0.9), 'g', 'filled');
        scatter(x_pie(prob_pie>0.9), -y_pie(prob_pie>0.9), 'm', 'filled');

        xlabel('Coordenada X');
        ylabel('Coordenada Y');
        title('Gráfico de los puntos');
        axis equal;
        xlim([0, max([x_manodch; x_manoizq; x_cola; x_pie])]);
        ylim([-max([y_manodch; y_manoizq; y_cola; y_pie]), 0]);

        legend(base_names.');
        grid on;

        hold off;

        nombre_grafico = fullfile(mypath, 'conducta','conducta.png');
        saveas(gcf, nombre_grafico);

        %quedarme solo con pvalor > 0.9
        % x_manodch(prob_manodch < 0.9) = NaN;
        % y_manodch(prob_manodch < 0.9) = NaN;
        % x_manoizq(prob_manoizq < 0.9) = NaN;
        % y_manoizq(prob_manoizq < 0.9) = NaN;
        % x_cola(prob_cola < 0.9) = NaN;
        % y_cola(prob_cola < 0.9) = NaN;
        % x_pie(prob_pie < 0.9) = NaN;
        % y_pie(prob_pie < 0.9) = NaN;

   
        % Calcula la derivada de las coordenadas en x e y
        dx_manodch = diff(x_manodch);
        dy_manodch = diff(y_manodch);

        % Grafica la derivada de las coordenadas
        figure;
        subplot(2, 1, 1);
        plot(dx_manodch,'b');
        xlabel('Tiempo');
        ylabel('Cambio en X');
        title(['Derivada de las coordenadas X de ', base_names{1}]);

        subplot(2, 1, 2);
        plot(dy_manodch,'b');
        xlabel('Tiempo');
        ylabel('Cambio en Y');
        title(['Derivada de las coordenadas Y de ', base_names{1}]);

        % Guarda la figura como un archivo de imagen PNG
        nombre_grafico = fullfile(mypath,'conducta', ['derivada',base_names{1},'.png']);
        saveas(gcf, nombre_grafico);

        % Calcula la derivada en ambas direcciones
        derivada_2d_dch = sqrt(dx_manodch.^2 + dy_manodch.^2);

        % Grafica la derivada en ambas direcciones
        figure;
        plot(derivada_2d_dch,'b');
        xlabel('Tiempo');
        ylabel('Cambio en posición (2D)');
        title(['Derivada de la posición de ', base_names{1},' 2D']);

        % Guarda la figura como un archivo de imagen PNG
        nombre_grafico = fullfile(mypath,'conducta', ['derivada_2D_',base_names{1},'.png']);
        saveas(gcf, nombre_grafico);

        % Calcula la derivada de las coordenadas en x e y
        dx_manoizq = diff(x_manoizq);
        dy_manoizq = diff(y_manoizq);

        % Grafica la derivada de las coordenadas
        figure;
        subplot(2, 1, 1);
        plot(dx_manoizq,'r');
        xlabel('Tiempo');
        ylabel('Cambio en X');
        title(['Derivada de las coordenadas X de ', base_names{2}]);

        subplot(2, 1, 2);
        plot(dy_manoizq,'r');
        xlabel('Tiempo');
        ylabel('Cambio en Y');
        title(['Derivada de las coordenadas Y de ', base_names{2}]);

        % Guarda la figura como un archivo de imagen PNG
        nombre_grafico = fullfile(mypath,'conducta', ['derivada',base_names{2},'.png']);
        saveas(gcf, nombre_grafico);

        % Calcula la derivada en ambas direcciones
        derivada_2d_izq = sqrt(dx_manoizq.^2 + dy_manoizq.^2);

        % Grafica la derivada en ambas direcciones
        figure;
        plot(derivada_2d_izq,'r');
        xlabel('Tiempo');
        ylabel('Cambio en posición (2D)');
        title(['Derivada de la posición de ', base_names{2},' 2D']);

        % Guarda la figura como un archivo de imagen PNG
        nombre_grafico = fullfile(mypath,'conducta', ['derivada_2D_',base_names{2},'.png']);
        saveas(gcf, nombre_grafico);

         % Calcula la derivada de las coordenadas en x e y
        dx_cola = diff(x_cola);
        dy_cola = diff(y_cola);

        % Grafica la derivada de las coordenadas
        figure;
        subplot(2, 1, 1);
        plot(dx_cola,'g');
        xlabel('Tiempo');
        ylabel('Cambio en X');
        title(['Derivada de las coordenadas X de ', base_names{3}]);

        subplot(2, 1, 2);
        plot(dy_cola,'g');
        xlabel('Tiempo');
        ylabel('Cambio en Y');
        title(['Derivada de las coordenadas Y de ', base_names{3}]);

        % Guarda la figura como un archivo de imagen PNG
        nombre_grafico = fullfile(mypath,'conducta', ['derivada',base_names{3},'.png']);
        saveas(gcf, nombre_grafico);

        % Calcula la derivada en ambas direcciones
        derivada_2d_cola = sqrt(dx_cola.^2 + dy_cola.^2);

        % Grafica la derivada en ambas direcciones
        figure;
        plot(derivada_2d_cola,'g');
        xlabel('Tiempo');
        ylabel('Cambio en posición (2D)');
        title(['Derivada de la posición de ', base_names{3},' 2D']);

        % Guarda la figura como un archivo de imagen PNG
        nombre_grafico = fullfile(mypath,'conducta', ['derivada_2D_',base_names{3},'.png']);
        saveas(gcf, nombre_grafico);

        % Calcula la derivada de las coordenadas en x e y
        dx_pie = diff(x_pie);
        dy_pie = diff(y_pie);

        % Grafica la derivada de las coordenadas
        figure;
        subplot(2, 1, 1);
        plot(dx_pie,'m');
        xlabel('Tiempo');
        ylabel('Cambio en X');
        title(['Derivada de las coordenadas X de ', base_names{4}]);

        subplot(2, 1, 2);
        plot(dy_pie,'m');
        xlabel('Tiempo');
        ylabel('Cambio en Y');
        title(['Derivada de las coordenadas Y de ', base_names{4}]);

        % Guarda la figura como un archivo de imagen PNG
        nombre_grafico = fullfile(mypath,'conducta', ['derivada',base_names{4},'.png']);
        saveas(gcf, nombre_grafico);

        % Calcula la derivada en ambas direcciones (norma euclidiana)
        derivada_2d_pie = sqrt(dx_pie.^2 + dy_pie.^2);

        % Grafica la derivada en ambas direcciones
        figure;
        plot(derivada_2d_pie,'m');
        xlabel('Tiempo');
        ylabel('Cambio en posición (2D)');
        title(['Derivada de la posición de ', base_names{4},' 2D']);

        % Guarda la figura como un archivo de imagen PNG
        nombre_grafico = fullfile(mypath,'conducta', ['derivada_2D_',base_names{4},'.png']);
        saveas(gcf, nombre_grafico);

        close all

        %4 gráficas y elección
        figure;
        subplot(4,1,1); %bodypart1
        plot(derivada_2d_dch,'b');
        xlabel('Tiempo');
        ylabel('Cambio en posición (2D)');
        title(base_names{1});
        subplot(4,1,2); %bodypart2
        plot(derivada_2d_izq,'r');
        xlabel('Tiempo');
        ylabel('Cambio en posición (2D)');
        title(base_names{2});
        subplot(4,1,3); %bodypart3
        plot(derivada_2d_cola,'g');
        xlabel('Tiempo');
        ylabel('Cambio en posición (2D)');
        title(base_names{3});
        subplot(4,1,4); %objectA
        plot(derivada_2d_pie,'m');
        xlabel('Tiempo');
        ylabel('Cambio en posición (2D)');
        title(base_names{4});

        nombre_grafico = fullfile(mypath,'conducta', 'derivadas_4partes.png');
        saveas(gcf, nombre_grafico);
   


        %derivada final
        figure;
        derivada_finalsuma = (derivada_2d_dch + derivada_2d_izq + derivada_2d_cola + derivada_2d_pie);
        plot(derivada_finalsuma)
        xlabel('Tiempo');
        ylabel('Derivada movimiento');
        title('Suma de movimiento de las 4 partes');
        grid on;

        % Guarda la figura como un archivo de imagen PNG en la carpeta "conducta"
        nombre_grafico = fullfile(mypath, 'conducta', 'suma_movimiento.png');
        saveas(gcf, nombre_grafico);

        % Calcula la diferencia entre las coordenadas de ambas manos
        dx_manos = x_manodch - x_manoizq;
        dy_manos = y_manodch - y_manoizq;

        % Calcula la derivada de la diferencia de las coordenadas de las manos en x e y
        ddx_manos = diff(dx_manos);
        ddy_manos = diff(dy_manos);

        % Calcula la norma euclidiana de las derivadas
        derivada_norma = sqrt(ddx_manos.^2 + ddy_manos.^2);

        % Grafica la derivada de las coordenadas
        figure;
        plot(derivada_norma,'Color', [0.35,0.13,0.39]);
        xlabel('Tiempo');
        ylabel('Derivada');
        title(['Derivada del movimiento de ',base_names{1},'-',base_names{2}]);
        legend('Derivada 2D');
        grid on;

        % Guarda la figura como un archivo de imagen PNG en la carpeta "conducta"
        nombre_grafico = fullfile(mypath, 'conducta', ['derivada_movimiento_',base_names{1},'-',base_names{2},'.png']);
        saveas(gcf, nombre_grafico);

        %Guardar todas las variables
        DLC.data.datos = datos;
        DLC.data.object1.dx = dx_manodch;
        DLC.data.object1.dy = dy_manodch;
        DLC.data.object1.twoD = derivada_2d_dch;
        DLC.data.object2.dx = dx_manoizq;
        DLC.data.object2.dy = dy_manoizq;
        DLC.data.object2.twoD = derivada_2d_izq;
        DLC.data.object3.dx = dx_cola;
        DLC.data.object3.dy = dy_cola;
        DLC.data.object3.twoD = derivada_2d_cola;
        DLC.data.object4.dx = dx_pie;
        DLC.data.object4.dy = dy_pie;
        DLC.data.object4.twoD = derivada_2d_pie;
        DLC.data.all = derivada_finalsuma;
        DLC.data.manos = derivada_norma;
        DLC.names = names;

    else
        disp('No se encontró ningún archivo CSV o se encontraron múltiples archivos CSV en la carpeta.');
        DLC.data.datos = NaN;
        DLC.data.both = NaN;
    end
end

% function DLC = openDLC(mypath)
%    % Obtiene la lista de archivos en la carpeta
%     archivos = dir(fullfile(mypath, '*.csv'));
% 
%     % Crea una carpeta para guardar las gráficas si no existe
%     if ~exist(fullfile(mypath, 'conducta'), 'dir')
%         mkdir(fullfile(mypath, 'conducta'));
%     end
% 
%     % Verifica si se encontró exactamente un archivo CSV
%     if numel(archivos) == 1
%         % Lee el archivo CSV como una tabla y asigna los títulos de las columnas
%         datos = readtable(fullfile(mypath, archivos(1).name), 'ReadVariableNames', false);
% 
%         % Interactuar con el usuario para ingresar nombres de columnas
%         prompt = {'bodypart 1', 'bodypart 2', 'bodypart 3', 'object A'};
%         dlg_title = 'Ingrese los nombres base para las categorías';
%         num_lines = 1;
%         defaultans = {'Mano dch', 'Mano izq', 'Cola','Pie'};
%         base_names = inputdlg(prompt,dlg_title,num_lines,defaultans);
% 
%         % Generar nombres de columnas basados en las bases ingresadas
%         column_names = cell(1, numel(base_names) * 3);
%         for i = 1:numel(base_names)
%             column_names{3*(i-1) + 1} = ['X ' base_names{i}];
%             column_names{3*(i-1) + 2} = ['Y ' base_names{i}];
%             column_names{3*(i-1) + 3} = ['prob ' base_names{i}];
%         end 
% 
%         titulos = column_names;
% % % Obtiene la lista de archivos en la carpeta
%     % archivos = dir(fullfile(mypath, '*.csv'));
%     % 
%     % % Crea una carpeta para guardar las gráficas si no existe
%     %     if ~exist(fullfile(mypath, 'conducta'), 'dir')
%     %         mkdir(fullfile(mypath, 'conducta'));
%     %     end
%     % 
%     % % Verifica si se encontró exactamente un archivo CSV
%     % if numel(archivos) == 1
%     %     % Define los títulos de las columnas
%     %     titulos = {'Frames','X mano dch','Y mano dch','prob mano dch','X mano izq','Y mano izq','prob mano izq','X cola','Y cola','prob cola','X pie','Y pie','prob pie'};
%     % 
%         % Lee el archivo CSV como una tabla y asigna los títulos a las columnas
%         datos = readtable(fullfile(mypath, archivos(1).name), 'ReadVariableNames', false);
%         datos.Properties.VariableNames = titulos;
% 
%         % Extrae las coordenadas de los puntos para cada objeto
%         x_manodch = datos{:, 2};
%         y_manodch = datos{:, 3};
%         x_manoizq = datos{:, 5};
%         y_manoizq = datos{:, 6};
%         x_cola = datos{:, 8};
%         y_cola = datos{:, 9};
%         x_pie = datos{:, 11};
%         y_pie = datos{:, 12};
% 
%         % Extrae las probabilidades de los puntos para cada objeto
%         prob_manodch = datos{:, 4};
%         prob_manoizq = datos{:, 7};
%         prob_cola = datos{:, 10};
%         prob_pie = datos{:, 13};
% 
%         % Grafica los puntos
%         figure;
%         hold on;
% 
%         scatter(x_manodch(prob_manodch>0.9), -y_manodch(prob_manodch>0.9), 'b', 'filled');
%         scatter(x_manoizq(prob_manoizq>0.9), -y_manoizq(prob_manoizq>0.9), 'r', 'filled');
%         scatter(x_cola(prob_cola>0.9), -y_cola(prob_cola>0.9), 'g', 'filled');
%         scatter(x_pie(prob_pie>0.9), -y_pie(prob_pie>0.9), 'm', 'filled');
% 
%         xlabel('Coordenada X');
%         ylabel('Coordenada Y');
%         title('Grafico de los puntos');
%         axis equal;
%         xlim([0, max([x_manodch; x_manoizq; x_cola; x_pie])]);
%         ylim([-max([y_manodch; y_manoizq; y_cola; y_pie]), -500]);
% 
%         legend('Mano derecha', 'Mano izquierda', 'Cola', 'Pie');
%         grid on;
% 
%         hold off;
% 
%         nombre_grafico = fullfile(mypath, 'conducta','conducta.png');
%         saveas(gcf, nombre_grafico);
% 
%         % Calcula la derivada de las coordenadas en x e y
%         dx_manodch = diff(x_manodch);
%         dy_manodch = diff(y_manodch);
% 
%         % Grafica la derivada de las coordenadas
%         figure;
%         subplot(2, 1, 1);
%         plot(dx_manodch);
%         xlabel('Tiempo');
%         ylabel('Cambio en X');
%         title('Derivada de las coordenadas X de la mano derecha');
% 
%         subplot(2, 1, 2);
%         plot(dy_manodch);
%         xlabel('Tiempo');
%         ylabel('Cambio en Y');
%         title('Derivada de las coordenadas Y de la mano derecha');
% 
%         % Guarda la figura como un archivo de imagen PNG
%         nombre_grafico = fullfile(mypath,'conducta', 'derivadasdcha.png');
%         saveas(gcf, nombre_grafico);
% 
%         % Calcula la derivada en ambas direcciones
%         derivada_2d_dch = sqrt(dx_manodch.^2 + dy_manodch.^2);
% 
%         % Grafica la derivada en ambas direcciones
%         figure;
%         plot(derivada_2d_dch);
%         xlabel('Tiempo');
%         ylabel('Cambio en posición (2D)');
%         title('Derivada de la posición de la mano derecha (2D)');
% 
%         % Guarda la figura como un archivo de imagen PNG
%         nombre_grafico = fullfile(mypath,'conducta', 'derivada_2D_dcha.png');
%         saveas(gcf, nombre_grafico);
% 
%         % Calcula la derivada de las coordenadas en x e y
%         dx_manoizq = diff(x_manoizq);
%         dy_manoizq = diff(y_manoizq);
% 
%         % Grafica la derivada de las coordenadas
%         figure;
%         subplot(2, 1, 1);
%         plot(dx_manoizq,'r');
%         xlabel('Tiempo');
%         ylabel('Cambio en X');
%         title('Derivada de las coordenadas X de la mano derecha');
% 
%         subplot(2, 1, 2);
%         plot(dy_manoizq,'r');
%         xlabel('Tiempo');
%         ylabel('Cambio en Y');
%         title('Derivada de las coordenadas Y de la mano derecha');
% 
%         % Guarda la figura como un archivo de imagen PNG
%         nombre_grafico = fullfile(mypath,'conducta', 'derivadasizq.png');
%         saveas(gcf, nombre_grafico);
% 
%         % Calcula la derivada en ambas direcciones
%         derivada_2d_izq = sqrt(dx_manoizq.^2 + dy_manoizq.^2);
% 
%         % Grafica la derivada en ambas direcciones
%         figure;
%         plot(derivada_2d_izq,'r');
%         xlabel('Tiempo');
%         ylabel('Cambio en posición (2D)');
%         title('Derivada de la posición de la mano derecha (2D)');
% 
%         % Guarda la figura como un archivo de imagen PNG
%         nombre_grafico = fullfile(mypath,'conducta', 'derivada_2D_izq.png');
%         saveas(gcf, nombre_grafico);
% 
%         % Calcula la diferencia entre las coordenadas de ambas manos
%         dx_manos = x_manodch - x_manoizq;
%         dy_manos = y_manodch - y_manoizq;
% 
%         % Calcula la derivada de la diferencia de las coordenadas de las manos en x e y
%         ddx_manos = diff(dx_manos);
%         ddy_manos = diff(dy_manos);
% 
%         % Calcula la norma euclidiana de las derivadas
%         derivada_norma = sqrt(ddx_manos.^2 + ddy_manos.^2);
% 
%         % Grafica la derivada de las coordenadas
%         figure;
%         plot(derivada_norma,'Color', [0.75, 0.1, 0.9]);
%         xlabel('Tiempo');
%         ylabel('Derivada');
%         title('Derivada del movimiento de ambas manos');
%         legend('Derivada 2D');
%         grid on;
% 
%         % Guarda la figura como un archivo de imagen PNG en la carpeta "conducta"
%         nombre_grafico = fullfile(mypath, 'conducta', 'derivada_movimiento_manos.png');
%         saveas(gcf, nombre_grafico);
% 
%     else
%         error('No se encontró ningún archivo CSV o se encontraron múltiples archivos CSV en la carpeta.');
%     end
% 
%     %Guardar todas las variables
%     DLC.data.datos = datos;
%     DLC.data.manodch.dx = dx_manodch;
%     DLC.data.manodch.dy = dy_manodch;
%     DLC.data.manodch.twoD = derivada_2d_dch;
%     DLC.data.manoizq.dx = dx_manoizq;
%     DLC.data.manoizq.dy = dy_manoizq;
%     DLC.data.manoizq.twoD = derivada_2d_izq;
%     DLC.data.derivManos = derivada_norma;

%     % Calcula la diferencia entre las coordenadas de ambas manos
    %     dx_manos = x_manodch - x_manoizq;
    %     dy_manos = y_manodch - y_manoizq;
    % 
    %     % Calcula la derivada de la diferencia de las coordenadas de las manos en x e y
    %     ddx_manos = diff(dx_manos);
    %     ddy_manos = diff(dy_manos);
    % 
    %     % Calcula la norma euclidiana de las derivadas
    %     derivada_norma = sqrt(ddx_manos.^2 + ddy_manos.^2);
    % 
    %     % Grafica la derivada de las coordenadas
    %     figure;
    %     plot(derivada_norma,'Color', [0.35,0.13,0.39]);
    %     xlabel('Tiempo');
    %     ylabel('Derivada');
    %     title(['Derivada del movimiento de ',base_names{1},'-',base_names{2}]);
    %     legend('Derivada 2D');
    %     grid on;
    % 
    %     % Guarda la figura como un archivo de imagen PNG en la carpeta "conducta"
    %     nombre_grafico = fullfile(mypath, 'conducta', ['derivada_movimiento_',base_names{1},'-',base_names{2},'.png']);
    %     saveas(gcf, nombre_grafico);
    % 
    %      % Calcula la norma euclidiana de las derivadas
    %     derivada_final = diff(derivada_norma);
    %     derivada_acel = sqrt(derivada_final.^2);
    % 
    %     % Grafica la derivada de las coordenadas
    %     figure;
    %     plot(derivada_acel,'Color', [0.35,0.13,0.39]);
    %     xlabel('Tiempo');
    %     ylabel('Derivada');
    %     title(['Aceleración del movimiento de ',base_names{1},'-',base_names{2}]);
    %     legend('Aceleración');
    %     grid on;
    % 
    % 
    %     % Guarda la figura como un archivo de imagen PNG en la carpeta "conducta"
    %     nombre_grafico = fullfile(mypath, 'conducta', ['aceleración_movimiento_',base_names{1},'-',base_names{2},'.png']);
    %     saveas(gcf, nombre_grafico);
    % 
    %     %Guardar todas las variables
    %     DLC.data.datos = datos;
    %     DLC.data.object1.dx = dx_manodch;
    %     DLC.data.object1.dy = dy_manodch;
    %     DLC.data.object1.twoD = derivada_2d_dch;
    %     DLC.data.object2.dx = dx_manoizq;
    %     DLC.data.object2.dy = dy_manoizq;
    %     DLC.data.object2.twoD = derivada_2d_izq;
    %     DLC.data.both = derivada_norma;
    %     DLC.data.acel = derivada_acel;
    %     DLC.names = names;
    % 
% 
% end

