function NAG()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NAG - Neural Analysis of Gain %
% Jon Boley                     %
% jdboley(a)purdue.edu          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function NAG(analysis)
% options for analysis:
%   'longrate'  (optimal gain restores log-term rate to normal)
%   'shortrate' (optimal gain restores short-term rate to normal)
%   'env'  (optimal gain restores env to normal)
%   'tfs'  (optimal gain restores tfs to normal)
% if nargin < 1, analysis='shortrate'; end;

GEN_SPIKES = 0;
ANALYZE_SPIKES = 0;

%% Initialization Routine
% Get Input Audio File Name   
PathName = 'C:\Audio Research\Audio\TIMIT\test\dr1\faks0\';
FileName = 'sa1.wav';
if (true) %set to false if you want to ask the user for an input file
    phonemeindx = textread([PathName FileName(1:end-3) 'phn'],'%*d %d %*s');
else
    [FileName,PathName] = uigetfile('*.wav',...
        'Select Input Audio File',[PathName FileName]);

    % Get TIMIT phoneme structure if needed
    reply = 1;%menu('Load TIMIT phoneme structure?','YES','NO');
    if (reply==1), % load index for end of each phoneme
        phonemeindx = textread([PathName FileName(1:end-3) 'phn'],'%*d %d %*s'); 
    end
    pause(0.1); % just wait for window to diappear
end

% Check for Neural Metrics Files
if strcmp(which('CCCanal'),'') % if we cannot find this function
    pathtool; % Let user add the appropriate directory to the path
    errordlg('Cannot Find Neural Metrics Files!');
end

%% Default Neural Model Params
lowCF=250; %Hz
highCF=4000; %Hz
numCFs=30;
spread=log2(highCF/lowCF); %number of octaves to span
CF_kHz=(250*2.^(0:spread/(numCFs-1):spread))/1000;
Cohc=1.0;  % initial OHC health
Cihc=1.0;  % initial IHC health
Nreps=120; % number of repetitions
ANmodel_Fs_Hz=100e3; % 100kHz sampling

%% Model hearing loss
FREQUENCIES = [250,500,1000,2000,3000,4000,6000]; % Hz
dBLoss = [20 20 30 40 45 50 50]; % Mild Hearing Loss (according to Bruce (ISAAR 2007)
impairment = 'mixed';
switch impairment
    case 'mixed'
        Dsd_OHC_Loss = 2/3*dBLoss;%zeros(1,length(FREQUENCIES));%24 *ones(1,length(FREQUENCIES));
    case 'ohc'
        Dsd_OHC_Loss = dBLoss;
    otherwise %case ihc
        Dsd_OHC_Loss = 0*dBLoss;
end
[Cohc_impaired,Cihc_impaired]=fitaudiogram(FREQUENCIES,dBLoss,Dsd_OHC_Loss);

% Display Model Information
file = dir('NAG.m');
disp(sprintf('Neural Analysis of Gain\nby Jon Boley'));
disp(sprintf('jdboley(a)purdue.edu\n%s',file.date));


%% Analyze the data
[input_orig,Fs] = readnist([PathName FileName]);
input_orig = input_orig/max(input_orig);

phones = 12;%1:length(phonemeindx);
levels = 70;%10:10:100;
gains = 0:5:5;%-40:5:40;
sponts = [50, 5, 0.25];
N_win_short = round(.000256*ANmodel_Fs_Hz); % 256us
N_win_long = round(.008*ANmodel_Fs_Hz); % 8ms
W_win_short = hann(N_win_short)/sum(hann(N_win_short));
W_win_long = hann(N_win_long)/sum(hann(N_win_long));
w_spont = [.6 .2 .2]; %[60%HSR 20%MSR 20%LSR]

%% initialization
% SynOutA = cell(length(gains),numCFs,length(sponts));
% SynOutB = cell(length(gains),numCFs,length(sponts));
% if (GEN_SPIKES)
%     SpikesA_plus = cell(length(gains),numCFs,length(sponts));
%     SpikesA_minus = cell(length(gains),numCFs,length(sponts));
%     SpikesB_plus = cell(length(gains),numCFs,length(sponts));
%     SpikesB_minus = cell(length(gains),numCFs,length(sponts));
% end
% if (ANALYZE_SPIKES)
%     Rate   = NaN*ones(length(gains),numCFs,length(sponts));
%     Difcor = NaN*ones(length(gains),numCFs,length(sponts));
%     Sumcor = NaN*ones(length(gains),numCFs,length(sponts));
%     Env = NaN*ones(length(gains),numCFs,length(sponts));
%     Tfs = NaN*ones(length(gains),numCFs,length(sponts));
%     OptimalGain = NaN*ones(length(levels));
% end


phone_index=1;
for phone=phones
    disp(sprintf('Processing phoneme #%d',phone));
    [input_model_phone,dBSPL_after,dur_sec] = InitPhoneme(input_orig,phonemeindx,phone,Fs,ANmodel_Fs_Hz);
    
    level_index=1;
    for OALevel_dBSPL=levels
        disp(sprintf('Processing level of %ddBSPL',OALevel_dBSPL));
        gain_index=1;        
        input_model_normal=InitLevel(input_model_phone,OALevel_dBSPL,dBSPL_after,ANmodel_Fs_Hz,dur_sec);

        for Gain_Adjust=gains
            disp(sprintf('Adjusting prescribed gain by %ddB',Gain_Adjust));
            
            for Fiber_Number=1:numCFs
                fprintf('.'); if (mod(Fiber_Number,5)==0), fprintf(' '); end
                if (Fiber_Number==numCFs), fprintf('\n'); end
                input_model = InitFiber(FREQUENCIES,dBLoss,CF_kHz,Fiber_Number,...
                    OALevel_dBSPL,Gain_Adjust,dBSPL_after,input_model_phone,ANmodel_Fs_Hz,dur_sec);

                spont_index=1;
                for SR_sps = sponts
                    %disp(sprintf('Processing spont rate of %d',SR_sps));
                    
                    [neurogramA1,neurogramA2,neurogramB1,neurogramB2] = ...
                        RunFibers(input_model,input_model_normal,gain_index,...
                        Fiber_Number,spont_index,CF_kHz,ANmodel_Fs_Hz,dur_sec,...
                        Cohc,Cihc,Cohc_impaired,Cihc_impaired,SR_sps,Nreps,...
                        GEN_SPIKES,ANALYZE_SPIKES,numCFs,w_spont);
                    
                    spont_index=spont_index+1;    
                end % end SR_sps
                
                % Filter neurogram data (one fiber at a time)
                neurogramA2(:,Fiber_Number) = filter(W_win_long,1,neurogramA1(:,Fiber_Number));
                neurogramB2(:,Fiber_Number) = filter(W_win_long,1,neurogramB1(:,Fiber_Number)); % apply long window
                neurogramA1(:,Fiber_Number) = filter(W_win_short,1,neurogramA1(:,Fiber_Number));
                neurogramB1(:,Fiber_Number) = filter(W_win_short,1,neurogramB1(:,Fiber_Number)); % apply short window
            end % end Fiber_Number
            
            % Save data to *.mat file
            note = datestr(now,'mmm_dd_yy');  % note to append to end of file name
            SaveData(phone,OALevel_dBSPL,impairment,Gain_Adjust,note);            
            
            gain_index=gain_index+1;
        end % end Gain_Adjust
        
        % reset metrics for next level calculation
        if (ANALYZE_SPIKES)
            clear Rate Difcor Sumcor Env Tfs
            Rate   = NaN*ones(length(gains),numCFs,length(sponts));
            Difcor = NaN*ones(length(gains),numCFs,length(sponts));
            Sumcor = NaN*ones(length(gains),numCFs,length(sponts));
            Env = NaN*ones(length(gains),numCFs,length(sponts));
            Tfs = NaN*ones(length(gains),numCFs,length(sponts));
        end
        
        level_index=level_index+1;
    end % end OALevel_dBSPL
    phone_index = phone_index+1;
end % end phone




% -------------------------
function [input_model_phone,dBSPL_after,dur_sec] = InitPhoneme(input_orig,phonemeindx,phone,Fs,ANmodel_Fs_Hz)
%%%%%%% Will need to include all previous phonemes
%     input = input_orig(1:phonemeindx(phone));

    % Just do one phoneme for now
    input=input_orig(phonemeindx(phone-1):phonemeindx(phone));
    
    % Resample for AN model
    dBSPL_before=20*log10(sqrt(mean(input.^2))/(20e-6));
    sfreq=Fs;
    sfreqNEW=ANmodel_Fs_Hz;
    P=round(sfreqNEW/10); Q=round(sfreq/10);  %Integers used to up sample
    if(P/Q*sfreq~=sfreqNEW), disp('Integer sfreq conversion NOT exact'); end
    Nfir=30;  % proportional to FIR filter length used for resampling: higher Nfir, better accuracy & longer comp time
    input_model_phone=resample(input,P,Q,Nfir);
    dBSPL_after=20*log10(sqrt(mean(input_model_phone.^2))/(20e-6));
    if abs(dBSPL_before-dBSPL_after)>2;
        error('RESAMPLING CHANGED input dBSPL by %f dB',dBSPL_after-dBSPL_before)
    end
    dur_sec = length(input_model_phone)/Fs;
%     clear input;


% -------------------------
function [neurogramA1,neurogramA2,neurogramB1,neurogramB2] = RunFibers(input_model,...
    input_model_normal,gain_index,Fiber_Number,spont_index,CF_kHz,ANmodel_Fs_Hz,dur_sec,Cohc,Cihc,...
    Cohc_impaired,Cihc_impaired,SR_sps,Nreps,GEN_SPIKES,ANALYZE_SPIKES,numCFs,w_spont,...
    neurogramA1,neurogramA2,neurogramB1,neurogramB2)

%% Run fiber (A+)
[timeout,meout,c1filterout,c2filterout,c1vihc,c2vihc,vihc,sout,psth500k] ...
    = zbcatmodel(input_model.',CF_kHz(Fiber_Number)*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.05,Cohc_impaired,Cihc_impaired,SR_sps);
SynOutA{gain_index,Fiber_Number,spont_index} = sout;
if (GEN_SPIKES)
    [sptimes nspikes]= SGfast([1/ANmodel_Fs_Hz, Nreps],sout);
    % Format spikes into NEL spikes format then cell array
    NELspikes=ANmodelSTs2nel(sptimes,Nreps);
    SpikeTrainsA_plus=nelSTs2cell(NELspikes);
    SpikesA_plus{gain_index,Fiber_Number,spont_index} = SpikeTrainsA_plus;
end
if (ANALYZE_SPIKES)
    Rate(gain_index,Fiber_Number,spont_index) = mean(sout); % SAVE RATE INFO
end

%% Run fiber (A-)
[timeout,meout,c1filterout,c2filterout,c1vihc,c2vihc,vihc,sout,psth500k] ...
    = zbcatmodel(-input_model.',CF_kHz(Fiber_Number)*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.05,Cohc_impaired,Cihc_impaired,SR_sps);
if (GEN_SPIKES)
    [sptimes nspikes]= SGfast([1/ANmodel_Fs_Hz, Nreps],sout);
    % Format spikes into NEL spikes format then cell array
    NELspikes=ANmodelSTs2nel(sptimes,Nreps);
    SpikeTrainsA_minus=nelSTs2cell(NELspikes);
    SpikesA_minus{gain_index,Fiber_Number,spont_index} = SpikeTrainsA_minus;
end

%% Run fiber (B+) - always normal
[timeout,meout,c1filterout,c2filterout,c1vihc,c2vihc,vihc,sout,psth500k] ...
    = zbcatmodel(input_model_normal.',CF_kHz(Fiber_Number)*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.05,Cohc,Cihc,SR_sps);
SynOutB{gain_index,Fiber_Number,spont_index} = sout;
if (GEN_SPIKES)
    [sptimes nspikes]= SGfast([1/ANmodel_Fs_Hz, Nreps],sout);
    % Format spikes into NEL spikes format then cell array
    NELspikes=ANmodelSTs2nel(sptimes,Nreps);
    SpikeTrainsB_plus=nelSTs2cell(NELspikes);
    SpikesB_plus{gain_index,Fiber_Number,spont_index} = SpikeTrainsB_plus;
end
%% Run fiber (B-) - always normal
[timeout,meout,c1filterout,c2filterout,c1vihc,c2vihc,vihc,sout,psth500k] ...
    = zbcatmodel(-input_model_normal.',CF_kHz(Fiber_Number)*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.05,Cohc,Cihc,SR_sps);
if (GEN_SPIKES)
    [sptimes nspikes]= SGfast([1/ANmodel_Fs_Hz, Nreps],sout);
    % Format spikes into NEL spikes format then cell array
    NELspikes=ANmodelSTs2nel(sptimes,Nreps);
    SpikeTrainsB_minus=nelSTs2cell(NELspikes);
    SpikesB_minus{gain_index,Fiber_Number,spont_index} = SpikeTrainsB_minus;
end
%%
if (ANALYZE_SPIKES)
    % Calculate env & tfs coding
    % Spike Analysis
    % Organize variables for CCCanal
    SpikeTrains=cell(2); % {condition (1,2), polarity (plus,minus)}
    SpikeTrains={SpikeTrainsA_plus,SpikeTrainsA_minus;SpikeTrainsB_plus,SpikeTrainsB_minus};

    % specify params to be used
    clear paramsIN
    paramsIN.durA_msec=dur_sec*1000;
    paramsIN.durB_msec=dur_sec*1000;
    paramsIN.CF_A_Hz=CF_kHz(Fiber_Number)*1000;
    paramsIN.CF_B_Hz=CF_kHz(Fiber_Number)*1000;
    paramsIN.MAXspikes=3000;
    paramsIN.PSD_LHfreqs_Hz=[0 64; 0 50];  %additional freq ranges to compute CCCenv for

    % [SACSCCfunctions,SACSCCmetrics,paramsOUT] = CCCanal(SpikeTrains,paramsIN,0);
    [SACSCCfunctions,SACSCCmetrics,paramsOUT] = CCCanal_4(SpikeTrains,paramsIN,0);

    Difcor(gain_index,Fiber_Number,spont_index) = SACSCCmetrics.DCpeak_A;
    Sumcor(gain_index,Fiber_Number,spont_index) = ...
        SACSCCmetrics.SCpeaks_A(find(strcmp('IFFTraw',SACSCCmetrics.SCpeaks_legend)));
    Env(gain_index,Fiber_Number,spont_index) = ...
        SACSCCmetrics.CCCenvs(find(strcmp('0-300, subBIAS',SACSCCmetrics.CCCenvs_legend)));
    Tfs(gain_index,Fiber_Number,spont_index) = SACSCCmetrics.CCCtfs;
end %if

%%
if ~isempty(SynOutA{gain_index,Fiber_Number,spont_index})
    neurogramA1(:,Fiber_Number) = neurogramA1(:,Fiber_Number) + ...
        w_spont(spont_index)*SynOutA{gain_index,Fiber_Number,spont_index}';
    neurogramB1(:,Fiber_Number) = neurogramB1(:,Fiber_Number) + ...
        w_spont(spont_index)*SynOutB{gain_index,Fiber_Number,spont_index}';
end


% -------------------------
function input_model_normal = InitLevel(input_model_phone,OALevel_dBSPL,dBSPL_after,ANmodel_Fs_Hz,dur_sec)

input_model_normal = input_model_phone*10^((OALevel_dBSPL-dBSPL_after)/20);
% REFIT and WINDOWwavefile at ANmodel_Fs
% Repeat or truncate waveform to fit requested stimulus duration:
input_model_normal = refit_waveform(input_model_normal,ANmodel_Fs_Hz,dur_sec*1000);
% Window waveform using linear rise/fall:
input_model_normal = window_waveform(input_model_normal,ANmodel_Fs_Hz,dur_sec*1000);
        
% Run Glasberg & Moore Compression
% input=RunAid(input,Fs,FREQUENCIES,NAL_IG)';
% dur_sec = length(input)/Fs;

% Apply OptimalGain to previous phonemes, and adjust gain of this phoneme.
% (NOTE: Will need to add attack/release here)
% EndIndex = floor(P/Q*phonemeindx(1));
% input_model(1:EndIndex) = input_model(1:EndIndex)*...
%         10^((OALevel_dBSPL-dBSPL_after+OptimalGain(1,level_index))/20);
% for PhoneNum=2:(phone-1)
%     StartIndex = floor(P/Q*phonemeindx(PhoneNum-1));
%     EndIndex = floor(P/Q*phonemeindx(PhoneNum));
%     input_model(StartIndex:EndIndex) = input_model(StartIndex:EndIndex)*...
%         10^((OALevel_dBSPL-dBSPL_after+OptimalGain(PhoneNum,level_index))/20);
% end


% ------------------------------
function input_model = InitFiber(FREQUENCIES,dBLoss,CF_kHz,Fiber_Number,...
    OALevel_dBSPL,Gain_Adjust,dBSPL_after,input_model_phone,ANmodel_Fs_Hz,dur_sec)

%% NEED TO ADD IN HUMAN HRTF
% REUG = [0 0 0 0 0 0 0]; % fake it for now

% % DSL4.0 table (from Dillon's book, p.243)
% % dBHL=0:5:110 (as defined by rows)
% %           250 500 750 1000    1500    2000    3000    4000    6000Hz
% DSL_REAG = [0   2   3   3       5       12      16      14      8;
%             3   4   5   5       8       15      18      17      11;
%             5   6   7   8       10      17      20      19      14;
%             7   8   10  10      13      19      23      21      17;
%             9   11  12  13      15      22      25      24      20;
%             12  13  14  15      18      24      28      27      23;
%             14  15  17  18      20      27      30      29      26;
%             17  18  19  21      23      30      33      32      29;
%             20  20  22  24      26      33      36      35      32;
%             22  23  25  27      29      36      39      38      36;
%             25  26  28  30      32      39      42      41      39;
%             29  29  31  33      35      42      45      45      43;
%             32  32  34  36      38      46      48      48      46;
%             36  35  37  40      42      49      52      51      50;
%             39  38  40  43      45      52      55      55      54;
%             43  42  43  46      48      56      59      58      58;
%             47  45  47  50      52      59      62      62      61;
%             51  48  50  53      55      63      66      65      65;
%             55  52  54  57      59      66      69      69      69;
%             59  55  57  60      62      70      73      73      0;
%             62  59  61  64      66      73      76      76      0;
%             0   62  64  68      70      77      80      80      0;
%             0   66  68  71      73      80      83      84      0];
% DSL_REAG=DSL_REAG(:,[1,2,4,6,7,8,9]); % ignore 750Hz and 1500Hz columns
% DSL_IG=zeros(length(FREQUENCIES),1);
% for i=1:length(FREQUENCIES)
%     DSL_IG(i) = DSL_REAG(floor(min(dBLoss(i),110)/5)+1,i) - REUG(i); % Calculate insertion gain
% end

% NAL equations
% frequency shaping (250,500,1000,2000,3000,4000,6000Hz)
k_NAL = [-18 -8 1 -1 -2 -2 -2]; % dB 
H_3FA = sum(dBLoss(2:4)) / 3; % sum up loss at 500Hz,1kHz,2kHz
X = 0.15*H_3FA;
R = 0.31; % NAL-RP formula  (not quite half-gain rule)
NAL_IG = X + R.*dBLoss + k_NAL; % insertion gain
NAL_IG = max(0,NAL_IG); % no negative gain

% Interpolate prescription
NAL = interp1(FREQUENCIES,NAL_IG,CF_kHz(Fiber_Number)*1000);
% DSL = interp1(FREQUENCIES,DSL_IG,CF_kHz(Fiber_Number)*1000);

adjustment = max(OALevel_dBSPL+NAL+Gain_Adjust-dBSPL_after,-dBSPL_after); % dB
% adjustment = max(OALevel_dBSPL+DSL+Gain_Adjust-dBSPL_after,-dBSPL_after);
input_model = input_model_phone*10^((adjustment)/20);

% REFIT and WINDOWwavefile at ANmodel_Fs
% Repeat or truncate waveform to fit requested stimulus duration:
input_model = refit_waveform(input_model,ANmodel_Fs_Hz,dur_sec*1000);
% Window waveform using linear rise/fall:
input_model = window_waveform(input_model,ANmodel_Fs_Hz,dur_sec*1000);

% ----------------------------
function [gain1,gain2,err1,err2] = OptimalGain(levels,gains,directory,impairment,date)
% [gain1,gain2,err1,err2] =
% OptimalGain(levels,gains,directory,impairment,date
% example:
% OptimalGain(30:10:100,-40:10:40,'archive\phone12\','mixed','Nov_24_08')

level_index=1;
for level=levels
    index=1;
    for gain = gains
        load([directory num2str(level) 'dBSPL\' impairment '\' num2str(gain) 'dBgain_' date], '-regexp', '^neurogram');
        err1(level_index,index) = mean(mean(abs(neurogramB1-neurogramA1))); % short window
        err2(level_index,index) = mean(mean(abs(neurogramB2-neurogramA2))); % long window
%         figure(1), subplot(length(gains),2,(2*index)-1), imagesc(neurogramB1'/max(max(neurogramB1)),[0 1]); axis off; title('Normal');
%         subplot(length(gains),2,2*index), imagesc(neurogramA1'/max(max(neurogramB1)),[0 1]); axis off; title(sprintf('Gain = %ddB',gain));
%         figure(2), subplot(length(gains),2,(2*index)-1), imagesc(neurogramB2'/max(max(neurogramB2)),[0 1]); axis off; title('Normal');
%         subplot(length(gains),2,2*index), imagesc(neurogramA2'/max(max(neurogramB2)),[0 1]); axis off; title(sprintf('Gain = %ddB',gain));
        clear -regexp ^neurogram
        index = index+1;
    end
    [C,I] = min(err1(level_index,:));
    gain1(level_index) = gains(I);
    [C,I] = min(err2(level_index,:));
    gain2(level_index) = gains(I);
    level_index=level_index+1;
    
%     pause;
end
figure, subplot(1,2,2), plot(levels,gain1,'x:'); 
axis([0 120 gains(1) gains(end)]); 
title(sprintf('NAL Gain Adjustment\nSpike Timing Information')); 
xlabel('Input Level (dB SPL)'); ylabel('Gain Adjustment (dB)');
subplot(1,2,1), plot(levels,gain2,'x:'); 
axis([0 120 gains(1) gains(end)]); 
title(sprintf('NAL Gain Adjustment\nAverage Discharge Rate')); 
xlabel('Input Level (dB SPL)'); ylabel('Gain Adjustment (dB)');


% ----------------------------
function SaveData(phone,OALevel_dBSPL,impairment,Gain_Adjust,note)
% Save file
MatFileName = sprintf('archive\\phone%d\\%1.0fdBSPL\\%s\\%1.0fdBgain_%s',...
    phone,OALevel_dBSPL,impairment,Gain_Adjust,note);
if ~exist(sprintf('archive\\phone%d\\%1.0fdBSPL\\%s',phone,OALevel_dBSPL,impairment),'dir')
    mkdir(sprintf('archive\\phone%d\\%1.0fdBSPL\\%s',phone,OALevel_dBSPL,impairment));
end
if (ANALYZE_SPIKES)
    save(MatFileName,'PathName','FileName','Difcor','Sumcor','Env',...
        'Tfs','OptimalGain','neurogramA1','neurogramA2');
else
    save(MatFileName,'PathName','FileName','neurogramA1','neurogramA2');
end
if Gain_Adjust==0, 
    save MatFileName neurogramB1 neurogramB2 -append
end

