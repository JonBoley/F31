% Calculate cross-CF correlations
% This script just calculates correlation re a specific CF
% (e.g., a CF centered on a vowel formant)

% The first fiber is the one centered on the feature of interest
SpikeTrains = cell(2); %{condition(1,2), polarity(plus,minus)}
SpikeTrains={Spikes_plus{1,LevelIndex,SNRindex},...
    Spikes_minus{1,LevelIndex,SNRindex};...
    Spikes_plus{FiberNumber,LevelIndex,SNRindex},...
    Spikes_minus{FiberNumber,LevelIndex,SNRindex}};

clear ParamsIN;
ParamsIN.durA_msec = dur_sec*1000;
ParamsIN.durB_msec = dur_sec*1000;
ParamsIN.CF_A_Hz = CF_kHz(1)*1000;
ParamsIN.CF_B_Hz = CF_kHz(FiberNumber)*1000;
ParamsIN.MAXspikes = 3000;
ParamsIN.PSD_LHfreqs_Hz = [0 64;0 50]; % additional freq ranges to calculate CCCenv for

[SACSCCfunctions,SACSCCmetrics,paramsOUT] = CCCanal_4(SpikeTrains,ParamsIN,0);

Rho(FiberNumber,LevelIndex,SNRindex) = ...
    SACSCCmetrics.SCCpeak_AB/sqrt(SACSCCmetrics.SACpeak_A*SACSCCmetrics.SACpeak_B);
CD(FiberNumber,LevelIndex,SNRindex) = ...
    SACSCCfunctions.delays_usec(SACSCCfunctions.SCC_AB_avg==max(SACSCCfunctions.SCC_AB_avg));

% % Manually determine CD (Correl,Delays,F0per_us,NSCC_CDs_usec,BFs_kHz,SCCind)
% [CD(FiberNumber,LevelIndex,SNRindex),NSCC_peaks{FiberNumber}] =...
% 	calcCD_manual(  SACSCCfunctions.SCC_AB_avg,...
%                     SACSCCfunctions.delays_usec,...
%                     1e6/formants(featureNum),...
%                     num2cell(CD(:,LevelIndex,SNRindex)),...
%                     mat2cell([CF_kHz(1) CF_kHz(FiberNumber)]),...
%                     1);

% Rho_env(FiberNumber) = SACSCCmetrics.CCCenv1;
% Rho_tfs(FiberNumber) = SACSCCmetrics.CCCtfs;


