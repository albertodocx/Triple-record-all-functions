function potenciasLFP(TDTdata, results,mypath)

pyr = results.lfp.datosfiltrados.';
pyrfs = TDTdata.streams.Wav1.fs;

% Análisis potencia:

pyrdata.dataLFP = double(pyr);
pyrdata.fsLFP = pyrfs;

pyrch_power = LFP_compute_psProfile(pyrdata);

% Extract specific data points from pyrch_power
    iugh = 'Theta';
    dataPoint2 = pyrch_power.infoRhs.theta.pmean  ;  
    dataPoint3 = pyrch_power.infoRhs.theta.pstd  ; 
    
    % Convert character vector to string
    dataPoint1 = string(iugh);

    % Create a table
    variableNames = {'Curva','Mean','Std'};
    dataTable = table(dataPoint1, dataPoint2, dataPoint3, 'VariableNames', variableNames);
    
    % Añadir otra fila con datos parecidos
    iugh2 = 'Delta';
    dataPoint2_2 = pyrch_power.infoRhs.delta.pmean;
    dataPoint3_2 = pyrch_power.infoRhs.delta.pstd;

    % Convert character vector to string
    dataPoint1_2 = string(iugh2);

    % Create the second row of the table
    newRow = table(dataPoint1_2, dataPoint2_2, dataPoint3_2, 'VariableNames', variableNames);

    % Añadir la segunda fila a la tabla existente
    dataTable = [dataTable; newRow];
    
    % Añadir otra fila con datos parecidos
    iugh3 = 'Gamma';
    dataPoint2_3 = pyrch_power.infoRhs.gamma.pmean;
    dataPoint3_3 = pyrch_power.infoRhs.gamma.pstd;

    % Convert character vector to string
    dataPoint1_3 = string(iugh3);

    % Create the second row of the table
    newRow_2 = table(dataPoint1_3, dataPoint2_3, dataPoint3_3, 'VariableNames', variableNames);

    % Añadir la segunda fila a la tabla existente
    dataTable = [dataTable; newRow_2];
    
    % Añadir otra fila con datos parecidos
    iugh4 = 'hfo';
    dataPoint2_4 = pyrch_power.infoRhs.hfo.pmean;
    dataPoint3_4 = pyrch_power.infoRhs.hfo.pstd;

    % Convert character vector to string
    dataPoint1_4 = string(iugh4);

    % Create the second row of the table
    newRow_3 = table(dataPoint1_4, dataPoint2_4, dataPoint3_4, 'VariableNames', variableNames);

    % Añadir la segunda fila a la tabla existente
    dataTable = [dataTable; newRow_3];
    
    % Añadir otra fila con datos parecidos
    iugh5 = 'Alpha';
    dataPoint2_5 = pyrch_power.infoRhs.alpha.pmean;
    dataPoint3_5 = pyrch_power.infoRhs.alpha.pstd;

    % Convert character vector to string
    dataPoint1_5 = string(iugh5);

    % Create the second row of the table
    newRow_4 = table(dataPoint1_5, dataPoint2_5, dataPoint3_5, 'VariableNames', variableNames);

    % Añadir la segunda fila a la tabla existente
    dataTable = [dataTable; newRow_4];

    % Display the updated table
    disp(dataTable);

    % Guardar la tabla en un archivo Excel en la carpeta seleccionada
    excelFileName = fullfile(mypath, 'Potencias.xlsx');
    writetable(dataTable, excelFileName);

    disp(['La tabla ha sido guardada en el archivo Excel: ' excelFileName]);
   

end
