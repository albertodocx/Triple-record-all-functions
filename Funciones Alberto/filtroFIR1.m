
function [signRip, signEnv] = filtroFIR1(Fs, matrix, pass1, pass2, order)
% INPUTS
% Fs: frecuencia de muestreo
% matrix: matriz a filtrar
% pass1: frecuencia mínima de paso
% pass2: frecuencia máxima de paso

% OUTPUTS
% signRip: filtrado de la matriz 'matrix' entre [pass1 pass2] 
% signEnv: envolvente del filtrado
if exist('order')
    N = order
else
    N  = 1700;  % Order
end
    Fstop1 = pass1-5;   % First Stopband 	
    Fpass1 = pass1;   % First Passband Frequency
    Fpass2 = pass2;  % Second Passband Frequency
    Fstop2 = pass2+5;  % Second Stopband Frequencyº
    Wstop1 = 1;     % First Stopband Weight
    Wpass  = 1;     % Passband Weight
    Wstop2 = 1;     % Second Stopband Weight
    if pass1 <15
        Fstop1= pass1-1;
    end
    if pass2 <15
        Fstop2= pass2+1;
    end
    % Calculate the coefficients using the FIRLS function. Ripple frequencies
    % filter.
    bRP  = firls(N, [0 Fstop1 Fpass1 Fpass2 Fstop2 Fs/2]/(Fs/2), [0 0 1 1 0 0], [Wstop1 Wpass Wstop2]);
    signRip = filtfilt(bRP,1,matrix);
    sign = 2*signRip.*signRip;
       
    % Amplitude of filtered signal (it's wide whenever there is a 70-250Hz event)
%     signEnv = movmean( movmean( sgolayfilt(sign,4,1001), pass1), pass2); 
%    signEnv =  sgolayfilt(sign,13,1001); 
    signEnv =  movmean( movmean( sgolayfilt(sign,11,1001), pass1), pass2);  % Me quedo con estos parámetros + movmean
    signEnv = abs(sqrt(signEnv));
    
    
end