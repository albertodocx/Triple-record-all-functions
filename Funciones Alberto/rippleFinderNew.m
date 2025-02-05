function [spont_rip, med, desv]=rippleFinderNew(data,ch, varargin)
% Inputs: 
% data: matrix with data
% ch: structure with channels
% chrip: channel for ripples
% Fs: Sampling frequency
% Threhold_ripple (Optional): threhold for ripple detection. Default 5
% Interval (Optional): In seconds. Vector of start and end of the interval
% apply nanstd (for detecting threshold). Default []

% Outputs:
% spont_rip: structure for spontaneous ripples
% med: mean of envelope pyr channel
% desv: std of envelope pyr channel

%% OPCIONES
% chrip,Fs,Fpass,treshold_ripple,interval
if any(strcmp(varargin,'Fpass'))
    Fpass = varargin{find(strcmp(varargin,'Fpass'))+1};
else
    Fpass = [90, 250];
end

if any(strcmp(varargin,'SD'))
    SD = varargin{find(strcmp(varargin,'SD'))+1};
end

if any(strcmp(varargin,'interval'))
    interval = varargin{find(strcmp(varargin,'interval'))+1};
else
    interval = [];
end

if any(strcmp(varargin,'threshold'))
    treshold_ripple = varargin{find(strcmp(varargin,'threshold'))+1};
else
    treshold_ripple = 5;
end

 if any(strcmp(varargin,'chrip'))
     chrip = varargin{find(strcmp(varargin,'chrip'))+1};    
 else
     chrip = ch.CA1;
 end

if any(strcmp(varargin,'Fs'))
    Fs = varargin{find(strcmp(varargin,'Fs'))+1};    
else
    Fs = 20000;
end

%% Select channels for the light pulse, PYR layer and RAD
% if ~isempty(ch.opto)
%     light=data(:,ch.opto);   
% else
%     light=zeros(size(data,1),1);
% end

% if max(light)>1
%     minimo=find(light(1:10*Fs)==min(light(1:10*Fs)),1,'first');  light=light-light(minimo,:); light(1:Fs)=0; light(end-Fs:end)=0; light=light/max(light);
% else
%     light=zeros(size(light));
% end
% CA1=data(:, ch.CA1); 
pyr=data(:,chrip);

 if ~isempty(ch.noise)
     chNoise=data(:,ch.noise); % Canal de corteza. Cuando hay ruido se vera afectado. Cuando hay ripple no
 else
     chNoise= zeros(size(data,1),1);
 end
% Areas where light On =
% pyr(find(light>2)-1000)=mean(pyr);  pyr(find(light>2)+1000)=mean(pyr); 
% chNoise(find(light>2)-1000)=mean(chNoise);  chNoise(find(light>2)+1000)=mean(chNoise); 

if ~isfield(ch, 'pip')
    pip=[];
else
    pip=data(:, ch.pip);    
end

if ~isfield(ch, 'Fpip')
    Fpip=[];
else
    Fpip=data(:, ch.Fpip);
end

% if isempty(ch.CA3)
%     CA3=[];
% else
%     CA3=data(:, ch.CA3);
% end


%% Filter ripple channel 300-3000 Hz to remove peaks
fprintf('Applying filter %d-%d to detect ripples... \n',300, 3000);  
if Fs < 6000
    df = designfilt('bandpassiir','FilterOrder',10,'HalfPowerFrequency1',300,'HalfPowerFrequency2',Fs/3,'SampleRate',Fs);

else
    df = designfilt('bandpassiir','FilterOrder',10,'HalfPowerFrequency1',300,'HalfPowerFrequency2',3000,'SampleRate',Fs);
end

pyr_filt = filtfilt(df,pyr-chNoise);


param = polyfit(chNoise, pyr,1);
if ~isempty(ch.noise)
    Y_fit = param(1) .*chNoise + param(2); %x405_2
    pyr_filt = filtfilt(df,pyr - Y_fit); 
else
    pyr_filt = 0;
end
%dF (units mV) is not dFF
%% Filter ripple channel 100-250 Hz
      
    tpoint = length(pyr);    % length of the signal in points
    fprintf('Applying filter %d-%d to detect ripples... \n',Fpass(1), Fpass(2));    
    
    [signRip, signEnv] = filtroFIR1(Fs, pyr-chNoise-pyr_filt, Fpass(1), Fpass(2));
    
% signEnv2 removes the data when the light pulse was applied (-1000 samples and +1000 samples)    

 signEnv2=signEnv;  
%  signEnv2(find(light>.2)-1000)=0;  signEnv2(find(light>.2)+1000)=0;  signEnv2(find(light>.2))=0;   signEnv2(signEnv2==0)=[];
 signEnv3=signEnv;  
% signEnv3(find(light>.2)-1000)=0;  signEnv3(find(light>.2)+1000)=0; signEnv3(find(light>.2))=0;
%  light2=light;   light2(find(light<.2))=0;  light2(find(light>=.2))=1; 
    for i=1:2
    signEnv2 = signEnv2 - mean(signEnv2);    
    end
 tpoint2 = length(signEnv2);
     
% if there is an interval, threshold is calculated only in that interval.
% Else, threshold is calculated from signEnv2 (excluding samples affected by the light pulse)
if ~exist('SD')
    if isempty(interval)
        thresh = nanstd(signEnv2)*treshold_ripple;
        
    else
        thresh = nanstd(signEnv(interval(1).*Fs:interval(2).*Fs))*treshold_ripple;
        
    end
else
    thresh = SD * treshold_ripple; 
end
    
% Find ripple candidates by searching local maximum samples in signEnv bigger than the threshold
k=1;   m= floor(Fs/20);    n= floor(Fs/200);
for i=(m+1):(tpoint-m)
        if ((signEnv3(i)>max(signEnv3(i-m:i-1))) && ...
                (signEnv3(i)>max(signEnv3(i+1:i+m))) && ...
                (signEnv3(i-n)>thresh) && (signEnv3(i+n)>thresh));
            events(k,1)=i;
            standard(k,1)= signEnv3(i)/(thresh/treshold_ripple);
            k=k+1;
        end
end
   
%% Remove repetitive candidates and candidates during the first second and last 3 seconds (avoid errors)
      
    events = sort(events); % unique(events);  %para que no se repitan. Lo enlimino
    for i=1:length(events)-1
       if (events(i+1)-events(i))< 0.1*Fs; % estaban 1600
           events(i)=0; 
           standard(i)=0;
       end
    end
    
    for i=1:length(events)
       if events(i)<Fs || events(i)> (length(pyr)-3*Fs) % Si se producen en el primer o 3 ultimos segundos los elimino (evitar errores de tamaño de matriz) 
           events(i)=0; 
           standard(i)=0;
       end
    end
    events(events==0)=[];
    standard(standard==0)=[];

%% Remove ripples events with high fast component
%fprintf('Applying filter %d-%d to detect ripples... \n',600, 3000);   
%[signRipFast, signEnvFast] = filtroFIR1(Fs, pyr-chNoise-pyr_filt, 600, 3000);
%signEnvFast=signEnvFast - min(signEnvFast); signEnvFast= signEnvFast/max(signEnvFast(1:10*Fs)); % Esta puesto en los 10 primeros seg
%signEnvFast(find(light>.2)-1000)=0;  signEnvFast(find(light>.2)+1000)=0; signEnvFast(find(light>.2))=0;
%signEnvRip=signEnv3 - min(signEnv3); signEnvRip=signEnvRip/max(signEnvRip(1:10*Fs)); % esta puesto en los 10 primeros seg

%events(signEnvFast(events)./signEnvRip(events)>=1)=[];
    
%% Separate spontaneous (1) and evoked (2) events: rip only type 1
% Esto por si se colara alguno. Antes hacia la separacion espontaneo y
% evocado en un mismo script

typeripple=[];

% if ~isempty(light2)
%     for i=1: length(events)
%         if light2(events(i))==0
%             typeripple(i)=1;
%         elseif light2(events(i))==1
%             typeripple(i)=2;
%         end
%     end
% else
%     typeripple(1: length(events),1) = 1;
% end
% rip=events(typeripple==1);
% std_rip=standard(typeripple==1);    

%% Plot spont ripples. Parte de Andrea
close all
fprintf('%d spontaneous events detected \n ',length(rip));
isripple=ones(size(rip));

if ~isempty(CA3)
    isripple = plotRipCA1_CA3(CA1, CA3, light, rip, Fs, [0.05, 0.05], 'Spontaneous events without light?');
else
    isripple = find_ripples_plot_Attila(pyr, light, rip, Fs, [0.05, 0.05], 'Spontaneous events without light?'); 
end

%% remove ripples manually discarted

rip(isripple==0)=[];
std_rip(isripple==0)=[];

if isempty(rip)
    std_rip=[];
    rip=[];
end

if ~isempty(CA3)
    isripCA3 = plotRipCA1_CA3(CA1, CA3, light, rip, Fs, [0.05, 0.05], 'Ripples in CA3?');
else
    isripCA3 = [];
end

%% Alignement by the minimum of the raw
        coil = floor(Fs * 0.01);
        minRaw = [];
        for i=1:length(rip)
            region = pyr((rip(i)-coil):(rip(i)+coil));
            [locmin] = islocalmin(region);
            pos= find(locmin);
            [~, minRaw] = min(region(locmin));  minRaw = pos(minRaw);
            % [~,minRaw] = min(pyr((rip(i)-coil):(rip(i)+coil)));
            rip(i)=rip(i)+(minRaw(1)-coil)-1;
        end    
%% Media y desviacion
if ~exist('SD')
if isempty(interval)
        desv=nanstd(signEnv2);  
        med=nanmean(signEnv2);
else
        desv=nanstd(signEnv(interval(1).*Fs:interval(2).*Fs));
        med=nanmean(signEnv2);
        %med=nanmean(signEnv(interval(1).*Fs:interval(2).*Fs));
end        
else 
    desv = SD;
    med = nanmean(signEnv2);
end
%% Choosing time for spont ripples (Output)
    spont_rip.idx_ripples= [];
    spont_rip.mat_CA1= [];
    spont_rip.mat_CA3= [];
    spont_rip.mat_pip= [];
    spont_rip.mat_Fpip= [];
    spont_rip.mat_light= [];
    spont_rip.mat_time=[];
    mat_env=[];
    mat_pyr=[];
% Spontaneous (time +- 500 ms)
    spont_rip.time_ripples= rip;  % Vector with time of max env of each ripple
    spont_rip.std_rip=std_rip;
    spont_rip.ripCA1_CA3=isripCA3;
 for i=1:length(rip)    
    spont_rip.idx_ripples(i,1)= Fs/2+1;
    spont_rip.mat_light(:,i)= light(rip(i) - Fs/2 : rip(i) + Fs/2);
    spont_rip.mat_time(:,i)= rip(i) - Fs/2 : rip(i) + Fs/2;
    mat_env(:,i)=signEnv(rip(i) - Fs/2 : rip(i) + Fs/2); 
    mat_pyr(:,i)= pyr(rip(i) - Fs/2 : rip(i) + Fs/2);
    
    if ~isempty(CA1)
        spont_rip.mat_CA1(:,i)= CA1(rip(i) - Fs/2 : rip(i) + Fs/2);     
    end
    
    if ~isempty(CA3)
        spont_rip.mat_CA3(:,i)= CA3(rip(i) - Fs/2 : rip(i) + Fs/2);     
    end
    
    if ~isempty(pip)
        spont_rip.mat_pip(:,i)= pip(rip(i) - Fs/2 : rip(i) + Fs/2);     
    end
    
    if ~isempty(Fpip)
        spont_rip.mat_Fpip(:,i)= Fpip(rip(i) - Fs/2 : rip(i) + Fs/2);     
    end
    
 end

% Duration of Spontaneus ripple
ripples=mat_pyr;
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