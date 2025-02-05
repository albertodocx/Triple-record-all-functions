%potencias
meanTheta_all = [];
meanDelta_all = [];
meanGamma_all = [];
meanHfo_all = [];
meanAlpha_all = [];

for i = 1:length(subfolders)
    disp(['Procesando subcarpeta número ', num2str(i)]);
    currentSubfolder = fullfile(mainFolder, subfolders(i).name);
    mypathTDT = strcat(currentSubfolder, '/fpdata.mat');
    mypathRs = strcat(currentSubfolder, '/Results.mat');
    load(mypathTDT, "TDTdata");
    load(mypathRs, "results");

    if ~exist(fullfile(currentSubfolder, 'definitivas'), 'dir')
        mkdir(fullfile(currentSubfolder, 'definitivas'));
    end

    subsubfolder = fullfile(currentSubfolder,'/definitivas');

    % Extraer datos de potencia
    pyr = results.lfp.datosfiltrados.';
    pyrfs = TDTdata.streams.Wav1.fs;
    pyrdata.dataLFP = double(pyr);
    pyrdata.fsLFP = pyrfs;
    pyrch_power = LFP_compute_psProfile(pyrdata);
    close all

    % Calcular las medias de potencia
    meanTheta = mean(pyrch_power.infoRhs.theta.pmean, 'all');
    meanDelta = mean(pyrch_power.infoRhs.delta.pmean, 'all');
    meanGamma = mean(pyrch_power.infoRhs.gamma.pmean, 'all');
    meanHfo = mean(pyrch_power.infoRhs.hfo.pmean, 'all');
    meanAlpha = mean(pyrch_power.infoRhs.alpha.pmean, 'all');

    % Almacenar las medias en matrices
    meanTheta_all = [meanTheta_all, meanTheta];
    meanDelta_all = [meanDelta_all, meanDelta];
    meanGamma_all = [meanGamma_all, meanGamma];
    meanHfo_all = [meanHfo_all, meanHfo];
    meanAlpha_all = [meanAlpha_all, meanAlpha];
end

% Crear un gráfico de líneas para las medias de cada potencia
figure(1);
subplot(2,1,1)
x_values = 1:length(subfolders);
plot(x_values, meanTheta_all, 'o-', 'DisplayName', 'Theta');
hold on;
plot(x_values, meanHfo_all, 'o-', 'DisplayName', 'HFO');
hold off;

% Personalizar el gráfico
xlabel('Subfolder');
ylabel('Valor Promedio');
title('Valor Promedio de Potencias para Diferentes Subfolders');
legend('Location', 'best');
xlim([0.5, length(subfolders) + 0.5]);

grid on;


%modulacion

media_pre = [];
sem_pre = [];
media_post = [];
sem_post = [];
rat_numbers = {}; % Inicializar array para almacenar los números de ratón



for i = 1:length(subfolders)
    disp(['vamos por i = ', num2str(i)]);
    currentSubfolder = fullfile(mainFolder, subfolders(i).name);
    mypathRs = fullfile(currentSubfolder, 'Results.mat');
    load(mypathRs, 'results')

    [~, numero_raton, ~] = fileparts(currentSubfolder);
    idx_ultimo_guion = regexp(numero_raton, '_[^_]*$');
    numero_raton = numero_raton(idx_ultimo_guion+1:end);
    rat_numbers{i} = numero_raton; % Almacenar el número de ratón

    dataFP_LFP = results.FP.analysis.PETH.onset.allintFP;

    denspreproc = table2array(dataFP_LFP(:, 3:4));
    discardindex = sum(isnan(denspreproc),2) == 0;
    denspreproc = denspreproc(discardindex,:);

    denspreFP_LFP = denspreproc(:, 1);
    denspostFP_LFP = denspreproc(:, 2);

    % Almacenar todos los datos pre y post
    all_denspreFP_LFP{i} = denspreFP_LFP;
    all_denspostFP_LFP{i} = denspostFP_LFP;

    % Calcular la media y el error estándar de la media (SEM)
    media_pre(end+1) = mean(denspreFP_LFP,'omitnan');
    sem_pre(end+1) = std(denspreFP_LFP) / sqrt(length(denspreFP_LFP));
    media_post(end+1) = mean(denspostFP_LFP,'omitnan');
    sem_post(end+1) = std(denspostFP_LFP) / sqrt(length(denspostFP_LFP));
end

% Crear un gráfico de dispersión para la media de pre y post
subplot(2,1,2)
hold on
for i = 1:length(media_pre)
    scatter(2*i-1, media_pre(i), 100, 'b', 'filled'); % Pre
    scatter(2*i, media_post(i), 100, 'r', 'filled');  % Post

    % Agregar barras de error (SEM)
    errorbar(2*i-1, media_pre(i), sem_pre(i), 'b', 'LineStyle', 'none', 'LineWidth', 1.5);
    errorbar(2*i, media_post(i), sem_post(i), 'r', 'LineStyle', 'none', 'LineWidth', 1.5);
end

% Etiquetas y título
ylabel('Densidad');
title('Media y SEM (pre y post)');
xticks(1:2:2*length(subfolders)); % Posiciones del eje x
xticklabels(rat_numbers); % Usar los números de ratón como etiquetas

xlim([0, 2*length(subfolders)+1]); % Límites del eje x
legend({'Pre', 'Post'}, 'Location', 'northwest');
grid on;