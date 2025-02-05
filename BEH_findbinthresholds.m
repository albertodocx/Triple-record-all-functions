
function [minthres, maxthres] = BEH_findbinthresholds(contvar, behfs)

    % Contvar: va a ser una matriz con todos los datos (todos los registros,
    % todos los animales, etc).
    % Cada fila es la n.

    % Output de min y max es un vector con los valores min y max de esa
    % variable para todos los registros/animales/etc.

    %     contvar = smooth(contvar, round(behfs), 'moving');
    
    contvarsm = nan(size(contvar));
    for ii = 1:size(contvar, 1)
        contvarsm(ii, 1:sum(~isnan(contvar))) = smooth(~isnan(contvar(ii)), round(behfs), 'moving');
    end

    meansm = mean(contvarsm, 2); stdsm = std(contvarsm, 2);
    minbin = min(contvar, [], "omitnan"); 
    maxbin = max(contvar, [], "omitnan"); 
    %
    % figure(1)
    % plot(minbin)
    % hold on
    % plot(maxbin)
    % ylim([-10, max(maxbin)])
    % hold off

    % Threshold in z-score:

    if min(minbin) == 0
        minthres = 0;
    elseif sum(minbin > 0) == length(minbin)
        minthres = min(minbin);
    else
        minthres = mean(minbin)-1.96*std(minbin);
    end

    %     maxthres = mean(maxbin)+1.96*std(maxbin); % 95% probability
    maxthres = maxbin;

end