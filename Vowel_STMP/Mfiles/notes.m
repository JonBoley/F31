% Jon's notes - when analyzing new data, do this:

%% Setup
setup_Vowel_STMP;

% Normal: 041811, 051211, 062711, 072011, 100611, 101111, 101711
% Impaired: 062311, 072111, 080111, 080911, 081511
% Aided: 050412
date='050412';

%% initial data processing
makeNOHRdataList(date);
storeBADlines(date);
% updateNOHR_DataList_excludeLines(date);

%% manually select SCC peaks
% batch_STMPplot_SCCfunctions;

%% Run cross-fiber analyses on all units
batch_run;
