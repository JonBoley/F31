% Jon's notes - when analyzing new data, do this:

%% Setup
setup_Vowel_STMP;

% Normal: 011811, 041811, 051211, 062711, 072011
% Impaired: 062311, 072111, 080111, 080911, 081511
% Aided: 
date='062711';

%% initial data processing
makeNOHRdataList(date);
storeBADlines(date);
% updateNOHR_DataList_excludeLines(date);

%% Run cross-fiber analyses on all units
batch_run;
