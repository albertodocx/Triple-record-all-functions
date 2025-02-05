function savefinal (mainFolder,subfolders,finalsubfolder)

[~, nombre_raton, ~] = fileparts(subfolders(1).name);
idx_medio_guion = regexp(nombre_raton, '(?<=_)[A-Za-z]\d+(?=_)', 'match', 'once');
nombre_raton = idx_medio_guion;

figure('Position',[218,47,4096,712]);
sgtitle(nombre_raton);

for i = 1:length(subfolders)
    disp(['vamos por i = ', num2str(i)]);
    currentSubfolder = fullfile(mainFolder, subfolders(i).name);
    mypathRs = strcat(currentSubfolder,'/Results.mat');
    load (mypathRs,"results")
    

    % Crea una carpeta para guardar las gráficas
    if ~exist(fullfile(currentSubfolder, 'definitivas'), 'dir')
        mkdir(fullfile(currentSubfolder, 'definitivas'));
    end

    subsubfolder = fullfile(currentSubfolder,'/definitivas');

    matrizfreqFP = (0:length(results.FP.Signals.DFFModZscore)-1)/results.FP.params.fs;
    matrizfreqLFP = (0:size(results.lfp.psprofilebueno.pS,2)-1)/2;
    matrizfreqcond = (0:length(results.DLC.data.manos.')-1)/15;

  % Extraer el nombre del ratón del nombre de la carpeta
    [~, numero_raton, ~] = fileparts(currentSubfolder);
    idx_ultimo_guion = regexp(numero_raton, '_[^_]*$');
    numero_raton = numero_raton(idx_ultimo_guion+1:end);



    %fotometría
    subplot (5,length(subfolders),(i));
    plot (matrizfreqFP,results.FP.Signals.DFFModZscore);
    xlim([0 matrizfreqFP(end)])
    title('Fotometría');
    % saveas(gcf,fullfile(subsubfolder,'1.jpg'));
    % pause (2);
    % close all
    hold on


    %electrofisiología
    fs = results.lfp.psprofilebueno.pSfreqs;
    T= matrizfreqLFP;
    PS=pow2db(results.lfp.psprofilebueno.pS);
    freqRange = (fs >= 0.5) & (fs <= 15);

    subplot (5,length(subfolders),(i+(length(subfolders))));
    imagesc(T, fs(freqRange), PS(freqRange,:));
    title('Espectrograma');
    colormap('jet');
    set(gca, 'YDir', 'normal');
    clim([0 100]);
    xlim([0 matrizfreqFP(end)])
    % saveas(gcf,fullfile(subsubfolder,'2.jpg'));
    % pause (2);
    % close all
    hold on


    %theta bouts
    boutsseg = results.lfp.thetaanalisis.bouts.segundos;
    subplot(5,length(subfolders),(i+(length(subfolders)*2)));
    for j = 1:size(boutsseg, 1)
        rectangle('Position', [boutsseg(j, 1), 0, boutsseg(j, 2) - boutsseg(j, 1), 1], 'FaceColor', 'k', 'EdgeColor', 'none');
        hold on;
    end

    title('Eventos de theta');
    xlim([0 matrizfreqFP(end)])
    ylim([0, 1]);
    grid on;
    % saveas(gcf,fullfile(subsubfolder,'3.jpg'));
    % pause(2)
    % close all
    hold on


    %conducta
    subplot(5,length(subfolders),(i+(length(subfolders)*3)));
    plot(matrizfreqcond,results.DLC.bouts.data.rescaled,"Color",[0.97,0.82,1.00]);
    hold on;
    plot(matrizfreqcond,results.DLC.bouts.data.smoothed,"Color",[0.35,0.13,0.39],"LineWidth",1.5);
    hold off;

    title('Conducta de movimiento');
    ylim([0, 1])
    xlim([0 matrizfreqFP(end)])
    hold on


    %mov bouts
    boutsdlc = results.DLC.bouts.segcut;
    subplot(5,length(subfolders),(i+(length(subfolders)*4)));
    for j = 1:size(boutsdlc, 1)
        rectangle('Position', [boutsdlc(j, 1), 0, boutsdlc(j, 2) - boutsdlc(j, 1), 1], 'FaceColor', 'k', 'EdgeColor', 'none');
        hold on;
    end

    title('Eventos de movimiento');
    xlabel('Tiempo(s)');
    ylim([0, 1]);
    xlim([0 matrizfreqFP(end)])
    text(mean(matrizfreqFP)-20,7.1,numero_raton,"Color",'red')
    grid on;
    hold on

end

saveas(gcf, fullfile(finalsubfolder, 'Global FP-LFP-Conducta.jpg'));
savefig(fullfile(finalsubfolder, 'Global FP-LFP-Conducta.fig'));
close all

for i = 1:length(subfolders)
    disp(['vamos por i = ', num2str(i)]);
    currentSubfolder = fullfile(mainFolder, subfolders(i).name);
    subsubfolder = fullfile(currentSubfolder,'/definitivas');
    mypathRs = strcat(currentSubfolder,'/Results.mat');
    load (mypathRs,"results")

    matrizfreqFP = (0:length(results.FP.Signals.DFFModZscore)-1)/results.FP.params.fs;
    matrizfreqLFP = (0:size(results.lfp.psprofilebueno.pS,2)-1)/2;
    matrizfreqcond = (0:length(results.DLC.data.manos.')-1)/15;

  % Extraer el nombre del ratón del nombre de la carpeta
    [~, numero_raton, ~] = fileparts(currentSubfolder);
    idx_ultimo_guion = regexp(numero_raton, '_[^_]*$');
    numero_raton = numero_raton(idx_ultimo_guion+1:end);

    figure('Position',[218,47,4096,712]);
    sgtitle(numero_raton);

    %fotometría
    subplot (5,1,1);
    plot (matrizfreqFP,results.FP.Signals.DFFModZscore);
    xlim([0 matrizfreqFP(end)])
    title('Fotometría');
    % saveas(gcf,fullfile(subsubfolder,'1.jpg'));
    % pause (2);
    % close all
    hold on


    %electrofisiología
    fs = results.lfp.psprofilebueno.pSfreqs;
    T= matrizfreqLFP;
    PS=pow2db(results.lfp.psprofilebueno.pS);
    freqRange = (fs >= 0.5) & (fs <= 15);

    subplot (5,1,2);
    imagesc(T, fs(freqRange), PS(freqRange,:));
    title('Espectrograma');
    colormap('jet');
    set(gca, 'YDir', 'normal');
    clim([0 100]);
    xlim([0 matrizfreqFP(end)])
    % saveas(gcf,fullfile(subsubfolder,'2.jpg'));
    % pause (2);
    % close all
    hold on


    %theta bouts
    boutsseg = results.lfp.thetaanalisis.bouts.segundos;
    subplot(5,1,3);
    for j = 1:size(boutsseg, 1)
        rectangle('Position', [boutsseg(j, 1), 0, boutsseg(j, 2) - boutsseg(j, 1), 1], 'FaceColor', 'k', 'EdgeColor', 'none');
        hold on;
    end

    title('Eventos de theta');
    xlim([0 matrizfreqFP(end)])
    ylim([0, 1]);
    grid on;
    % saveas(gcf,fullfile(subsubfolder,'3.jpg'));
    % pause(2)
    % close all
    hold on


    %conducta
    subplot(5,1,4);
    plot(matrizfreqcond,results.DLC.bouts.data.rescaled,"Color",[0.97,0.82,1.00]);
    hold on;
    plot(matrizfreqcond,results.DLC.bouts.data.smoothed,"Color",[0.35,0.13,0.39],"LineWidth",1.5);
    hold off;

    title('Conducta de movimiento');
    ylim([0, 1])
    xlim([0 matrizfreqFP(end)])
    hold on


    %mov bouts
    boutsdlc = results.DLC.bouts.segcut;
    subplot(5,1,5);
    for j = 1:size(boutsdlc, 1)
        rectangle('Position', [boutsdlc(j, 1), 0, boutsdlc(j, 2) - boutsdlc(j, 1), 1], 'FaceColor', 'k', 'EdgeColor', 'none');
        hold on;
    end

    title('Eventos de movimiento');
    xlabel('Tiempo(s)');
    ylim([0, 1]);
    xlim([0 matrizfreqFP(end)])
    grid on;
    hold off
    saveas(gcf, fullfile(subsubfolder, ['Global FP-LFP-Conducta',numero_raton,'.jpg']));
    savefig(fullfile(subsubfolder, ['Global FP-LFP-Conducta',numero_raton,'.fig']));
    close all

end