function currentsignals = conducta(TDTdata)

Fotoisos = TDTdata.streams.x405N.data;
Fotodata = TDTdata.streams.x465N.data;
Electro = TDTdata.streams.Wav1.data;
TTL = TDTdata.streams.Fie1.data;

% Establecer un umbral para detectar cruces de TTL
umbralCruce = 1;

% Encontrar cruces ascendentes y descendentes
crucesAscendentes = TTL > umbralCruce;
crucesDescendentes = TTL < umbralCruce;

% Identificar los índices donde se produce un cruce ascendente
indicesCrucesAscendentes = find(diff(crucesAscendentes) == 1) + 1;

% Identificar los índices donde se produce un cruce descendente
indicesCrucesDescendentes = find(diff(crucesDescendentes) == 1) + 1;

% Asegurarse de que hay al menos un cruce descendente después de cada cruce ascendente
indicesPicos = [];
for i = 1:length(indicesCrucesAscendentes)
    indiceAscendente = indicesCrucesAscendentes(i);
    indicesDescendentesDespues = indicesCrucesDescendentes(indicesCrucesDescendentes > indiceAscendente);

    if ~isempty(indicesDescendentesDespues)
        indiceDescendente = indicesDescendentesDespues(1);
        indicesPicos = [indicesPicos, round((indiceAscendente + indiceDescendente) / 2)];
    end
end

% Identificar el primer y último pico
primerPico = indicesPicos(1);
ultimoPico = indicesPicos(end);

% Asegurarse de que todas las señales tengan la misma longitud
Fotoisos = Fotoisos(primerPico:ultimoPico);
Fotodata = Fotodata(primerPico:ultimoPico);
Electro = Electro(primerPico:ultimoPico);
TTL = TTL(primerPico:ultimoPico);
Electro = Electro*1000;

% Actualizar la estructura TDTdata
currentsignals.Fotoisos = Fotoisos;
currentsignals.Fotodata = Fotodata;
currentsignals.Electro= Electro;
currentsignals.TTL = TTL;

% Crear un vector de tiempo para el eje x (asumiendo una frecuencia de muestreo de 1 Hz)
frecuenciaMuestreo = 1;  % ajusta según tu frecuencia de muestreo
tiempo = (0:length(TTL)-1) / frecuenciaMuestreo;

% Mostrar resultados
disp(['La señal TTL tiene ' num2str(length(indicesPicos)) ' frames.']);

% Plotear las señales después del corte
figure;
subplot(4, 1, 1);
plot(tiempo, Fotoisos, 'LineWidth', 2);
xlabel('Tiempo (s)');
ylabel('Fotoisos');
title('Señal Fotoisos después del Corte');
grid on;

subplot(4, 1, 2);
plot(tiempo, Fotodata, 'LineWidth', 2);
xlabel('Tiempo (s)');
ylabel('Fotodata');
title('Señal Fotodata después del Corte');
grid on;

subplot(4, 1, 3);
plot(tiempo, Electro, 'LineWidth', 2);
xlabel('Tiempo (s)');
ylabel('Electro');
title('Señal Electro después del Corte');
grid on;

subplot(4, 1, 4);
stairs(tiempo, TTL, 'LineWidth', 2);
hold on;
xlabel('Tiempo (s)');
ylabel('Valor de la Señal TTL');
title('Detección de Picos en la Señal TTL después del Corte');
grid on;
hold off;


end
