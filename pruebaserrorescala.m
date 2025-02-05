basicData = LFPdata;
Fs = fs;
nch = 16;
tWin = 2;
dtWin = 0.5; 

A = pS(:, :, 1);

sum(A(freq>=hfo_cf(1)&freq<=hfo_cf(2), win))

length(freq>=hfo_cf(1)&freq<=hfo_cf(2))

sum(freq>=hfo_cf(1)&freq<=hfo_cf(2))

which LCN_compute_psProfile.m


data =  LFPdata;

factor_ds = 10;
disp('downsampling...'); data_ds = downsample(data,factor_ds);
fs_ds = fs/factor_ds;

[pS,band,freq]=powerSpectrumProfileAlb(LFPdata(:,1:16),fs,1); 

data = results.LFPds.dsdata;
d = results.LFPds.dsdata;
Fs = 2000;
timeSpec = 1;
 

dprueba = d(:, 1);
media = mean(dprueba(1:10000));
tiledlayout(2, 1)
nexttile
plot(dprueba(1:10000), 'b')
nexttile
plot(dprueba(1:10000) - media, 'm')

datos = dprueba(1:10000) - media;
mean(datos)
mean(dprueba)

d2 = d - mean(d,1);

newloadLFPdata = LFPdata;
newresults = results;
original = results_rawLFP;


