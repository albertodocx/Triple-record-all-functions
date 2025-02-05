function saveanalisisFPLFP(mainFolder, subfolders)

[~, nombre_raton, ~] = fileparts(subfolders(1).name);
idx_medio_guion = regexp(nombre_raton, '(?<=_)[A-Za-z]\d+(?=_)', 'match', 'once');
nombre_raton = idx_medio_guion;

media_pre = [];
sem_pre = [];
media_post = [];
sem_post = [];
rat_numbers = {}; % Inicializar array para almacenar los números de ratón

% Inicializar arrays para almacenar todos los datos pre y post
all_denspreFP_LFP = {};
all_denspostFP_LFP = {};

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

% Crear una tabla con los datos pre
T_pre = table(rat_numbers', all_denspreFP_LFP', 'VariableNames', {'Rat_Number', 'Density_Pre'});

% Crear una tabla con los datos post
T_post = table(rat_numbers', all_denspostFP_LFP', 'VariableNames', {'Rat_Number', 'Density_Post'});

% Escribir las tablas de datos pre y post en archivos Excel
filename_pre = fullfile(mainFolder, ['AnalysisResults_Pre_',nombre_raton,'Fotometría-LFP.xlsx']);
writetable(T_pre, filename_pre);

filename_post = fullfile(mainFolder, ['AnalysisResults_Post_',nombre_raton,'Fotometría-LFP.xlsx']);
writetable(T_post, filename_post);

% Combinar los resultados en una sola tabla (opcional)
T_raw_data = [T_pre(:, 1:2), T_post(:, 2)];
filename_raw_data = fullfile(mainFolder, ['AnalysisResults_RawData_',nombre_raton,'Fotometría-LFP.xlsx']);
writetable(T_raw_data, filename_raw_data);

% Crear un gráfico de dispersión para la media de pre y post
figure('Position',[145,334,1003,420])
hold on
for i = 1:length(media_pre)
    scatter(2*i-1, media_pre(i), 100, 'MarkerFaceColor',[0.50,0.00,1.00],'MarkerEdgeColor','none'); % Pre
    scatter(2*i, media_post(i), 100, 'MarkerFaceColor',[0.99,0.50,0.03],'MarkerEdgeColor','none');  % Post

    % Agregar barras de error (SEM)
    errorbar(2*i-1, media_pre(i), sem_pre(i), 'Color',[0.50,0.00,1.00], 'LineStyle', 'none', 'LineWidth', 1.5);
    errorbar(2*i, media_post(i), sem_post(i), 'Color',[0.99,0.50,0.03], 'LineStyle', 'none', 'LineWidth', 1.5);
end

% Etiquetas y título
ylabel('Densidad');
title('Media y SEM (pre y post) - Fotometría//LFP');
xticks(1:2:2*length(subfolders)); % Posiciones del eje x
xticklabels(rat_numbers); % Usar los números de ratón como etiquetas

xlim([0, 2*length(subfolders)+1]); % Límites del eje x
legend({'Pre', 'Post'}, 'Location', 'northwest');
grid on;
saveas (gcf,fullfile(mainFolder,['Modulación ', nombre_raton,'-FotometríayLFP.jpg']));
end








% 
%     dataFP_LFP = results.FP.analysis.PETH.onset  .allintFP;
%     %dataFP_Conducta = results.DLC.conducta.PETH.onset.allintFP; 
%     %datatheta_LFP = results.lfp.lfp_zscore.ththbouts.PETH.onset.allintFP;
% 

