function [spont_rip, med, desv]= rippleFinder_1ch (data, varargin)

% Outputs:
% spont_rip: structure for spontaneous ripples
% med: mean of envelope pyr channel
% desv: std of envelope pyr channel

%% Default setting for the optional parameters   
p = inputParser;
addParameter(p,'Fs', 20000, @isnumeric)  % Fs: Sampling frequency
addParameter(p,'Fpass', [70, 250], @isnumeric)
addParameter(p,'Tresrip', 3, @isnumeric) % Threhold_ripple (Optional): threhold for ripple detection. Default 3
addParameter(p,'Treswave', 2, @isnumeric) 
addParameter(p,'Interval', [], @isnumeric) % Interval (Optional): In seconds. Vector of start and end of the interval
addParameter(p,'Signrip', [], @isnumeric) 
addParameter(p,'Signenv', [], @isnumeric) 
addParameter(p,'Signlow', [], @isnumeric) 

parse(p, varargin{:})
% Variables
Fs = p.Results.Fs;
Fpass = p.Results.Fpass;
interval = p.Results.Interval;
treshold_ripple= p.Results.Tresrip;
treshold_wave = p.Results.Treswave;
signRip = p.Results.Signrip;
signEnv = p.Results.Signenv;
signLow = p.Results.Signlow;

%% Filter ripple channel 100-250 Hz      
    tpoint = length(data);    % length of the signal in points
   
    if isempty(signRip) || isempty(signEnv)
        fprintf('Applying filter %d-%d to detect ripples... \n',Fpass(1), Fpass(2)); 
        [signRip, signEnv] = filtroFIR1(Fs, data, Fpass(1), Fpass(2));  
    end

    if isempty(signLow)
        signLow = filtroLOW(Fs, data, 8);   
    end

    signLow = abs(signLow);

% if there is an interval, threshold is calculated only in that interval.
% Else, threshold is calculated from signEnv2 (excluding samples affected by the light pulse)
    if isempty(interval)
        thresh = nanstd(signEnv)*treshold_ripple;    
        threshLOW = nanstd(signLow)*treshold_wave;
    else
        thresh = nanstd(signEnv(interval(1).*Fs:interval(2).*Fs))*treshold_ripple;    
        threshLOW = nanstd(signLow(interval(1).*Fs:interval(2).*Fs))*treshold_wave;
    end
    
% Find ripple candidates by searching local maximum samples in signEnv bigger than the threshold
k=1;   m= floor(Fs/20);    n= floor(Fs/200);
for i=(m+1):(tpoint-m)
        if ((signEnv(i)>max(signEnv(i-m:i-1))) && ...
                (signEnv(i)>max(signEnv(i+1:i+m))) && ...
                (signEnv(i-n)>thresh) && (signEnv(i+n)>thresh) && ...
                signLow(i-n)>threshLOW) && (signLow(i+n)>threshLOW);
            events(k,1)=i;
            standard(k,1)= signEnv(i)/(thresh/treshold_ripple);
            k=k+1;
        end
end
   
%% Remove repetitive candidates and candidates during the first second and last second (avoid errors)
      
    events = sort(events); % unique(events);  %para que no se repitan. Lo enlimino
    for i=1:length(events)-1
       if (events(i+1)-events(i))< 0.1*Fs; % estaban 1600
           events(i)=0; 
           standard(i)=0;
       end
    end
    
    for i=1:length(events)
       if events(i)<Fs || events(i)> (length(data)-Fs) % Si se producen en el primer ultimo segundo los elimino (evitar errores de tamaño de matriz) 
           events(i)=0; 
           standard(i)=0;
       end
    end
    events(events==0)=[];
    standard(standard==0)=[];
   
%% Separate spontaneous (1) and evoked (2) events: rip only type 1
% Esto por si se colara alguno. Antes hacia la separacion espontaneo y
% evocado en un mismo script
rip=events;
std_rip=standard;    

%% Plot spont ripples. Parte de Andrea
close all
fprintf('%d spontaneous events detected \n ',length(rip));
isripple=ones(size(rip));

isripple = find_ripples_plot(data, signLow, rip, Fs, [0.15, 0.15], 'Spontaneous events without light?'); 


%% remove ripples manually discarted

rip(isripple==0)=[];
std_rip(isripple==0)=[];

if isempty(rip)
    std_rip=[];
    rip=[];
end

%% Alignement by the minimum of the raw
        coil = floor(Fs * 0.01);
        minRaw = [];
        for i=1:length(rip)
            region = data((rip(i)-coil):(rip(i)+coil));
            [locmin] = islocalmin(region);
            pos= find(locmin);
            [~, minRaw] = min(region(locmin));  minRaw = pos(minRaw);
            % [~,minRaw] = min(pyr((rip(i)-coil):(rip(i)+coil)));
            rip(i)=rip(i)+(minRaw(1)-coil)-1;
        end    
%% Media y desviacion
if isempty(interval)
        desv=nanstd(signEnv);  
        med=nanmean(signEnv);
else
        desv=nanstd(signEnv(interval(1).*Fs:interval(2).*Fs));
        med=nanmean(signEnv);        
end        
        
%% Choosing time for spont ripples (Output)
    spont_rip.idx_ripples= [];
    spont_rip.mat_CA1= [];    
    spont_rip.mat_time=[];
    mat_env=[];    
% Spontaneous (time +- 500 ms)
    spont_rip.time_ripples= rip;  % Vector with time of max env of each ripple
    spont_rip.std_rip=std_rip;
    
 for i=1:length(rip)    
    spont_rip.idx_ripples(i,1)= Fs/2+1;    
    spont_rip.mat_time(:,i)= rip(i) - Fs/2 : rip(i) + Fs/2;
    mat_env(:,i)=signEnv(rip(i) - Fs/2 : rip(i) + Fs/2);     
    spont_rip.mat_CA1(:,i)= data(rip(i) - Fs/2 : rip(i) + Fs/2);            
 end

% Duration of Spontaneus ripple
ripples= spont_rip.mat_CA1;
centro=ceil(size(ripples,1)./2);
inicio=[]; final=[];
if isempty(ripples)
    spont_rip.duraRipple=[];
else        
    deriEnv =diff(diff(mat_env)); deriEnv=[zeros(2,size(deriEnv,2));deriEnv].*1000;
    
    for i=1: size(mat_env,2)
        inicio(i,1)=find(deriEnv(centro-0.075*Fs : centro-0.01*Fs,i)==max(deriEnv(centro-0.075*Fs : centro-0.01*Fs,i)),1,'first');
        final(i,1)=find(deriEnv(centro+0.01*Fs : centro+0.075*Fs,i)==max(deriEnv(centro+0.01*Fs : centro+0.075*Fs,i)),1,'first')+0.085*Fs;
    end
    
    spont_rip.duraRipple=(final-inicio)./(Fs./1000);
    
end

   
end