
%%%%%%%%%%%%%%%%%%%%%%   PIPELINE FOR LFP-ANALYSIS %%%%%%%%%%%%%%%%%%%%%

%% 1. Load and save all files

myfolder = uigetdir();
ICL_loadLFPfiles(myfolder)

% 1.1. Apply 50 Hz notch filter (to process gamma later):

ICL_notch50filter(myfolder)

%% 2. Select channels of interest:

ICL_idLFPch(myfolder)


%% 3. Detect and get measures for theta events:

% thetadata2 = ICL_thetaevent(myfolder);
thetadata2 = ICL_thetaevent_zscore(myfolder);

%% 4. Detect and get measures for ripple events (classic filter; CNN not implemented yet): 

rippledata = ICL_getripples(myfolder);

