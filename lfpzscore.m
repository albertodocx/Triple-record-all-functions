function lfp_zscore = lfpzscore(results, mypath)

% Crear la carpeta "lfp" si aún no existe
if ~exist(fullfile(mypath, 'lfp'), 'dir')
    mkdir(fullfile(mypath, 'lfp'));
end
data = results.lfp.psprofilebueno.pS;
fs = results.lfp.psprofilebueno.pSfreqs;
datadb = pow2db(data);

% % Zcore curvas de potencia

delta_cf = [1, 3]; % Delta power: between 1 and 3Hz
theta_cf = [6, 9]; % Theta power: between 4 and 12Hz
alpha_cf = [9, 16]; % Alpha power: between 9 and 16Hz
gamma_cf = [20, 90]; % Gamma power: between 20 and 90Hz
hfo_cf = [100, 500]; % HFO power: between 100 and 500Hz

% Theta
fstheta = (fs>=theta_cf(1)&fs<=theta_cf(2));
thetadata = datadb(fstheta,:);
vectortheta = mean(thetadata,1);
thetazscore = zscore(vectortheta);
%thetaenv = sgolayfilt(smooth(thetazscore),51,1729);
thetaenv = movmean( movmean( sgolayfilt(thetazscore,50,51), 20), 20);

meanTheta = mean (thetazscore);
stdTheta = std (thetazscore);
meanThetaenv = mean (thetaenv);
stdThetaenv = std (thetaenv);

% Guardar la figura de Theta en la carpeta "lfp"
figure('Position',[46,468,1323,286]);
plot (smooth(thetazscore), "Color",[0.3,0.72,1.00]);
hold on
plot (thetaenv,"Color",[0.00,0.45,0.74],"LineWidth",1.5)
title ('Theta Zscore')
nombre_grafico = fullfile(mypath, 'lfp', 'Theta Zscore.png');
saveas(gcf, nombre_grafico);
close all; % Cerrar la figura

lfp_zscore.data.thetazscore = thetazscore;
lfp_zscore.power.thetazscore.mean = meanTheta;
lfp_zscore.power.thetazscore.std = stdTheta;
lfp_zscore.power.thetazscore.meanenv = meanThetaenv;
lfp_zscore.power.thetazscore.stdenv = stdThetaenv;

% Delta
fsdelta = (fs>=delta_cf(1)&fs<=delta_cf(2));
deltadata = datadb(fsdelta,:);
vectordelta = mean(deltadata,1);
deltazscore = zscore(vectordelta);
%deltaenv = sgolayfilt(smooth(deltazscore),51,1729);
deltaenv = movmean( movmean( sgolayfilt(deltazscore,50,51), 20), 20);

meanDelta = mean (deltazscore);
stdDelta = std (deltazscore);

% Guardar la figura de Delta en la carpeta "lfp"
figure('Position',[46,468,1323,286]);
plot (smooth(deltazscore), "Color",[1.00,0.60,0.43]);
hold on
plot (deltaenv,"Color",[0.85,0.33,0.10],"LineWidth",1.5)
title ('Delta Zscore')
nombre_grafico = fullfile(mypath, 'lfp', 'Delta Zscore.png');
saveas(gcf, nombre_grafico);
close all; % Cerrar la figura

lfp_zscore.data.deltazscore = deltazscore;
lfp_zscore.power.deltazscore.mean = mean (deltazscore);
lfp_zscore.power.deltazscore.std = std (deltazscore);
lfp_zscore.power.deltazscore.meanenv = mean (deltaenv);
lfp_zscore.power.deltazscore.stdenv = std (deltaenv);

% Gamma
fsgamma = (fs>=gamma_cf(1)&fs<=gamma_cf(2));
gammadata = datadb(fsgamma,:);
vectorgamma = mean(gammadata,1);
gammazscore = zscore(vectorgamma);
%gammaenv = sgolayfilt(smooth(gammazscore),51,1729);
gammaenv = movmean( movmean( sgolayfilt(gammazscore,50,51), 20), 20);

meanGamma = mean (gammazscore);
stdGamma = std (gammazscore);

% Guardar la figura de Gamma en la carpeta "lfp"
figure('Position',[46,468,1323,286]);
plot (smooth(gammazscore), "Color",[0.98,0.83,0.49]);
hold on
plot (gammaenv,"Color",[0.93,0.69,0.13],"LineWidth",1.5)
title ('Gamma Zscore')
nombre_grafico = fullfile(mypath, 'lfp', 'Gamma Zscore.png');
saveas(gcf, nombre_grafico);
close all; % Cerrar la figura

lfp_zscore.data.gammazscore = gammazscore;
lfp_zscore.power.gammazscore.mean = mean (gammazscore);
lfp_zscore.power.gammazscore.std = std (gammazscore);
lfp_zscore.power.gammazscore.meanenv = mean (gammaenv);
lfp_zscore.power.gammazscore.stdenv = std (gammaenv);

% HFO
fshfo = (fs>=hfo_cf(1)&fs<=hfo_cf(2));
hfodata = datadb(fshfo,:);
vectorhfo = mean(hfodata,1);
hfozscore = zscore(vectorhfo);
%hfoenv = sgolayfilt(smooth(hfozscore),51,1729);
hfoenv = movmean( movmean( sgolayfilt(hfozscore,50,51), 20), 20);

meanHfo = mean (hfozscore);
stdHfo = std (hfozscore);

% Guardar la figura de HFO en la carpeta "lfp"
figure('Position',[46,468,1323,286]);
plot (smooth(hfozscore), "Color",[0.91,0.51,1.00]);
hold on
plot (hfoenv,"Color",[0.49,0.18,0.56],"LineWidth",1.5)
title ('HFO Zscore')
nombre_grafico = fullfile(mypath, 'lfp', 'Hfo Zscore.png');
saveas(gcf, nombre_grafico);
close all; % Cerrar la figura

lfp_zscore.data.hfozscore = hfozscore;
lfp_zscore.power.hfozscore.mean = mean (hfozscore);
lfp_zscore.power.hfozscore.std = std (hfozscore);
lfp_zscore.power.hfozscore.meanenv = mean (hfoenv);
lfp_zscore.power.hfozscore.stdenv = std (hfoenv);

% Alpha
fsalpha = (fs>=alpha_cf(1)&fs<=alpha_cf(2));
alphadata = datadb(fsalpha,:);
vectoralpha = mean(alphadata,1);
alphazscore = zscore(vectoralpha);
%alphaenv = sgolayfilt(smooth(alphazscore),51,1729);
alphaenv = movmean( movmean( sgolayfilt(alphazscore,50,51), 20), 20);

meanAlpha = mean (alphazscore);
stdAlpha = std (alphazscore);

% Guardar la figura de Alpha en la carpeta "lfp"
figure('Position',[46,468,1323,286]);
plot (smooth(alphazscore), "Color",[0.69,0.91,0.41]);
hold on
plot (alphaenv,"Color",[0.47,0.67,0.19],"LineWidth",1.5)
title ('Alpha Zscore')
nombre_grafico = fullfile(mypath, 'lfp', 'Alpha Zscore.png');
saveas(gcf, nombre_grafico);
close all; % Cerrar la figura

lfp_zscore.data.alphazscore = alphazscore;
lfp_zscore.power.alphazscore.mean = mean (alphazscore);
lfp_zscore.power.alphazscore.std = std (alphazscore);
lfp_zscore.power.alphazscore.meanenv = mean (alphaenv);
lfp_zscore.power.alphazscore.stdenv = std (alphaenv);

%Create a table
variableNames = {'Curva','Mean','Std'};
dataTable = table("Theta", meanTheta, stdTheta, 'VariableNames', variableNames);
newRow = table("Delta", meanDelta, stdDelta, 'VariableNames', variableNames);
newRow_2 = table("Gamma", meanGamma, stdGamma, 'VariableNames', variableNames);
newRow_3 = table("hfo", meanHfo, stdHfo, 'VariableNames', variableNames);
newRow_4 = table("Alpha", meanAlpha, stdAlpha, 'VariableNames', variableNames);
dataTable = [dataTable;newRow;newRow_2;newRow_3;newRow_4];
disp(dataTable);

% Plotting mean and standard deviation separately
figure;
plot(1:numel(meanTheta), meanTheta,"Marker",'o',"LineStyle",'-',"Color",[0.00,0.45,0.74],"MarkerFaceColor",[0.00,0.45,0.74],"MarkerSize",10);
hold on;
plot(1:numel(meanDelta), meanDelta,"Marker",'o',"LineStyle",'-',"Color",[0.85,0.33,0.10],"MarkerFaceColor",[0.85,0.33,0.10],"MarkerSize",10);
plot(1:numel(meanGamma), meanGamma,"Marker",'o',"LineStyle",'-',"Color",[0.93,0.69,0.13],"MarkerFaceColor",[0.93,0.69,0.13],"MarkerSize",10);
plot(1:numel(meanHfo), meanHfo,"Marker",'o',"LineStyle",'-',"Color",[0.49,0.18,0.56],"MarkerFaceColor",[0.49,0.18,0.56],"MarkerSize",10);
plot(1:numel(meanAlpha), meanAlpha,"Marker",'o',"LineStyle",'-',"Color",[0.47,0.67,0.19],"MarkerFaceColor",[0.47,0.67,0.19],"MarkerSize",10);
hold off

% % Personalizar el gráfico
xlabel('Curva');
ylabel('Valor Promedio');
title('Valor Promedio para Diferentes Curvas');
legend('Theta', 'Delta', 'Gamma', 'hfo', 'Alpha');
nombre_grafico = fullfile(mypath, 'lfp', 'Potencias medias.png');
saveas(gcf, nombre_grafico);
close all


% Coincidencia entre eventos y zscored theta
boutsseg = results.lfp.thetaanalisis.bouts.muestras;

% Calcular el rango de valores de thetaenv
min_theta = min(thetaenv);
max_theta = max(thetaenv);
max_range = max(abs(min_theta), abs(max_theta)); % Obtener el valor absoluto máximo

% Calcular los límites en y para centrar alrededor de 0
y_limit = [-max_range, max_range];

% Crear la figura
figure('Position',[46,468,1323,286]);


% Trazar el segundo gráfico primero
for i = 1:size(boutsseg, 1)
    rectangle('Position', [boutsseg(i, 1), y_limit(1), boutsseg(i, 2) - boutsseg(i, 1), diff(y_limit)], 'FaceColor', [0.8,0.8,0.8], 'EdgeColor', 'none');
    hold on;
end
xdddlim = [1, numel(thetaenv)];
line(xdddlim, [0, 0], 'Color', 'k', 'LineStyle', '--');

% Trazar el primer gráfico encima del segundo
plot(thetaenv, 'Color', [0.00,0.45,0.74],"LineWidth",1.5);
ylim(y_limit); % Establecer límites centrados en 0
xlim([1, numel(thetaenv)]); % Establecer límites en x
xlabel('Tiempo');
ylabel('Valor');
title('Gráfico combinado de Thetaenv y Eventos');
grid on;
nombre_grafico = fullfile(mypath, 'lfp', 'ThetaZscore_Eventos.png');
saveas(gcf, nombre_grafico);
pause (2)
close all

% Coincidencia de potencia de theta con los theta-bouts
DFF = thetazscore;
fsDFF = 2;
eventos = floor((results.lfp.thetaanalisis.bouts.segundos)*fsDFF);
fsev = fsDFF;
maxlength = floor(length(DFF)/fsDFF);
lfppath = fullfile(mypath, 'lfp');

lfp_zscore.ththbouts = BEH_PETHonsetv2(DFF, fsDFF, eventos, fsev, 'dffmovement', lfppath,'maxlength',maxlength,'AUCint',[-7 0; 0 7], 'Pre', 7, 'Post', 7, 'bin', 0.2);

end