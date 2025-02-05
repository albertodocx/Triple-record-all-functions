
function [minDist] = CFC_calcMDB(Location, fs)


% 5. Distancia mínima al borde (el que sea) en cada punto del registro.
% Output: vector de misma longitud que pos con el valor en cm de la
% distancia mínima al borde. 

rightlimx = 525; % Límites definitivos 
leftlimx = 125;
highlimy = 440;
lowlimy = 40;

% Location = BHD.Location;
locx = Location(:, 1);
locy = Location(:, 2);
locx(locx > rightlimx) = rightlimx;
locx(locx < leftlimx) = leftlimx;
locy(locy > highlimy) = highlimy;
locy(locy < lowlimy) = lowlimy;
Location = [locx locy];

rightdist = rightlimx-Location(:, 1);
leftdist = Location(:, 1)-leftlimx;
highdist = highlimy - Location(:, 2);
lowdist = Location(:, 2)-lowlimy;

alldist = [rightdist, leftdist,highdist,lowdist].';
frompxtocm = 16;
minDist = min(alldist)/frompxtocm;

maxd = 12.5; 
mind = 0; 
minDist(minDist < mind) = 0;
minDist(minDist > maxd) = maxd;

threeminidx = round(fs*60*3);
% minDist3 = minDist(1:threeminidx);


end