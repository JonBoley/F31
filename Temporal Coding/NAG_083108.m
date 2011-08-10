%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NAG - Neural Analysis of Gain %
% Jon Boley                     %
% jdboley(a)purdue.edu          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; home;

%% Initialization Routine

% Get Input Audio File Name    
[FileName,PathName,FilterIndex] = uigetfile('*.wav','Select Input Audio File','C:\Audio\TIMIT\test\dr1\faks0\sa1.wav');

% Get TIMIT phoneme structure if needed
reply = input('Load TIMIT phoneme structure? Y/N [Y]: ', 's');
if isempty(reply), reply = 'Y'; end
if strcmp(upper(reply),'Y'), % load index for end of each phoneme
    phonemeindx = textread([PathName FileName(1:end-3) 'phn'],'%*d %d %*s'); 
end

% Check for Neural Metrics Files
if strcmp(which('CCCanal'),'') % if we cannot find this function
    pathtool; % Let user add the appropriate directory to the path
    errordlg('Cannot Find Neural Metrics Files!');
end

% Neural Model Params
lowCF=250; %Hz
highCF=4000; %Hz
numCFs=20;
spread=log2(highCF/lowCF); %number of octaves to span
CF_kHz=(250*2.^(0:spread/(numCFs-1):spread))/1000;
SR_sps=50; % spontaneous rate
Cohc=1.0;  % initial OHC health
Cihc=1.0;  % initial IHC health
Nreps=120; % number of repetitions
ANmodel_Fs_Hz=100e3; % 100kHz sampling

% Model hearing loss
FREQUENCIES = [250,500,1000,2000,3000,4000,6000]; % Hz
dBLoss = [20 20 30 40 45 50 50]; % Mild Hearing Loss (according to Bruce (ISAAR 2007)
Dsd_OHC_Loss = 24 *ones(1,length(FREQUENCIES));
[Cohc_impaired,Cihc_impaired,OHC_Loss]=fitaudiogram(FREQUENCIES,dBLoss,Dsd_OHC_Loss);

% NEED TO ADD IN HUMAN HRTF
REUG = [0 0 0 0 0 0 0]; % fake it for now

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

% Display Model Information
file = dir('NAG.m');
disp(sprintf('Neural Analysis of Gain\nby Jon Boley'));
disp(sprintf('jdboley(a)purdue.edu\n%s',file.date));


%% Analyze the data
[input_orig,Fs] = readnist([PathName FileName]);
input_orig = input_orig/max(input_orig);

phones = 12;%1:length(phonemeindx);
levels = 100:-10:30;%10:10:100;
gains = 0;%-40:5:40;
sponts = 50;%[0.25, 5, 50];

% initialization
Rate   = NaN*ones(length(phones),length(levels),length(gains),numCFs,length(sponts));
Difcor = NaN*ones(length(phones),length(levels),length(gains),numCFs,length(sponts));
Sumcor = NaN*ones(length(phones),length(levels),length(gains),numCFs,length(sponts));
Env = NaN*ones(length(phones),length(levels),length(gains),numCFs,length(sponts));
Tfs = NaN*ones(length(phones),length(levels),length(gains),numCFs,length(sponts));
Rate_Normal   = NaN*ones(length(phones),length(levels),numCFs,length(sponts));
Difcor_Normal = NaN*ones(length(phones),length(levels),numCFs,length(sponts));
Sumcor_Normal = NaN*ones(length(phones),length(levels),numCFs,length(sponts));
Env_Normal = NaN*ones(length(phones),length(levels),numCFs,length(sponts));
Tfs_Normal = NaN*ones(length(phones),length(levels),numCFs,length(sponts));
phone_index=1;
for phone=phones
    disp(sprintf('Processing phoneme #%d',phone));
%     if phone==1,
%         input = input_orig(1:phonemeindx(phone));
%     else
%         input = input_orig(phonemeindx(phone-1):phonemeindx(phone));
%     end
    input=input_orig;
    dur_sec = length(input)/Fs;
    
    
    level_index=1;
    for OALevel_dBSPL=levels
        disp(sprintf('Processing level of %ddBSPL',OALevel_dBSPL));
        gain_index=1;
        
        input=RunAid(input,FREQUENCIES,NAL_IG);
        
        for Gain_Adjust=gains
            disp(sprintf('Adjusting prescribed gain by %ddB',Gain_Adjust));
            for Fiber_Number=1:numCFs
                disp(sprintf('Processing fiber #%d of %d',Fiber_Number,numCFs));
                spont_index=1;
                for SR_sps = sponts
                    disp(sprintf('Processing spont rate of %d',SR_sps));
                    % Interpolate prescription
                    NAL = interp1(FREQUENCIES,NAL_IG,CF_kHz(Fiber_Number)*1000);
        %             DSL = interp1(FREQUENCIES,DSL_IG,CF_kHz(Fiber_Number)*1000);

                    % Resample for AN model
                    dBSPL_before=20*log10(sqrt(mean(input.^2))/(20e-6));
                    sfreq=Fs;
                    sfreqNEW=ANmodel_Fs_Hz;
                    P=round(sfreqNEW/10); Q=round(sfreq/10);  %Integers used to up sample
                    if(P/Q*sfreq~=sfreqNEW), disp('Integer sfreq conversion NOT exact'); end
                    Nfir=30;  % proportional to FIR filter length used for resampling: higher Nfir, better accuracy & longer comp time
                    input_model=resample(input,P,Q,Nfir);
                    dBSPL_after=20*log10(sqrt(mean(input_model.^2))/(20e-6));
                    if abs(dBSPL_before-dBSPL_after)>2;
                        error(sprintf('RESAMPLING CHANGED input dBSPL by %f dB',dBSPL_after-dBSPL_before))
                    end
%                     clear input;

                    numreps=1;
                    if 0;%Gain_Adjust==gains(1) % if we're on the first set of calculations
                        numreps=2;  % also run normal
                    end
        
                    for i=1:numreps
                        adjustment=0;
                        % Apply gain (prescription + adjustment)
                        if i==1 % aided+impaired
                            adjustment = max(OALevel_dBSPL+(dBLoss(1)/2)-dBSPL_after,-dBSPL_after); % dB
%                             adjustment = max(OALevel_dBSPL+NAL+Gain_Adjust-dBSPL_after,-dBSPL_after); % dB
%                             adjustment = max(OALevel_dBSPL+DSL+Gain_Adjust-dBSPL_after,-dBSPL_after);
                        elseif i==2 % normal
                            adjustment = max(OALevel_dBSPL-dBSPL_after,-dBSPL_after);
                        end
                        input_model = input_model*10^((adjustment)/20);
                        
                        % REFIT and WINDOWwavefile at ANmodel_Fs
                        % Repeat or truncate waveform to fit requested stimulus duration:
                        input_model = refit_waveform(input_model,ANmodel_Fs_Hz,dur_sec*1000);
                        % Window waveform using linear rise/fall:
                        input_model = window_waveform(input_model,ANmodel_Fs_Hz,dur_sec*1000);         

                        % Run fiber (A+)
                        if i==1 % aided+impaired
                            [timeout,meout,c1filterout,c2filterout,c1vihc,c2vihc,vihc,sout,psth500k] ...
                                = zbcatmodel(input_model.',CF_kHz(Fiber_Number)*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.05,Cohc_impaired,Cihc_impaired,SR_sps);
                        else % normal
                            [timeout,meout,c1filterout,c2filterout,c1vihc,c2vihc,vihc,sout,psth500k] ...
                                = zbcatmodel(input_model.',CF_kHz(Fiber_Number)*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.05,Cohc,Cihc,SR_sps);
                        end
                        [sptimes nspikes]= SGfast([1/ANmodel_Fs_Hz, Nreps],sout);
                        % Format spikes into NEL spikes format then cell array
                        NELspikes=ANmodelSTs2nel(sptimes,Nreps);
                        SpikeTrainsA_plus=nelSTs2cell(NELspikes);

                        if i==1 % aided+impaired
                            Rate(phone_index,level_index,gain_index,Fiber_Number,spont_index) = mean(sout); % SAVE RATE INFO
                        elseif i==2 % normal
                            Rate_Normal(phone_index,level_index,Fiber_Number,spont_index) = mean(sout);
                        end

                        % Run fiber (A-)
                        if i==1 % aided+impaired
                            [timeout,meout,c1filterout,c2filterout,c1vihc,c2vihc,vihc,sout,psth500k] ...
                                = zbcatmodel(-input_model.',CF_kHz(Fiber_Number)*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.05,Cohc_impaired,Cihc_impaired,SR_sps);
                        else % normal
                            [timeout,meout,c1filterout,c2filterout,c1vihc,c2vihc,vihc,sout,psth500k] ...
                                = zbcatmodel(-input_model.',CF_kHz(Fiber_Number)*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.05,Cohc,Cihc,SR_sps);
                        end
                        [sptimes nspikes]= SGfast([1/ANmodel_Fs_Hz, Nreps],sout);
                        % Format spikes into NEL spikes format then cell array
                        NELspikes=ANmodelSTs2nel(sptimes,Nreps);
                        SpikeTrainsA_minus=nelSTs2cell(NELspikes);
                        
                        % Run fiber (B+) - always normal
                        [timeout,meout,c1filterout,c2filterout,c1vihc,c2vihc,vihc,sout,psth500k] ...
                            = zbcatmodel(input_model.',CF_kHz(Fiber_Number)*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.05,Cohc,Cihc,SR_sps);  
                        [sptimes nspikes]= SGfast([1/ANmodel_Fs_Hz, Nreps],sout);
                        % Format spikes into NEL spikes format then cell array
                        NELspikes=ANmodelSTs2nel(sptimes,Nreps);
                        SpikeTrainsB_plus=nelSTs2cell(NELspikes);
                        % Run fiber (B-) - always normal
                        [timeout,meout,c1filterout,c2filterout,c1vihc,c2vihc,vihc,sout,psth500k] ...
                            = zbcatmodel(-input_model.',CF_kHz(Fiber_Number)*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.05,Cohc,Cihc,SR_sps);
                        [sptimes nspikes]= SGfast([1/ANmodel_Fs_Hz, Nreps],sout);
                        % Format spikes into NEL spikes format then cell array
                        NELspikes=ANmodelSTs2nel(sptimes,Nreps);
                        SpikeTrainsB_minus=nelSTs2cell(NELspikes);

                        % Calculate env & tfs coding
                        %% Spike Analysis
                        % Organize variables for CCCanal
                        SpikeTrains=cell(2); % {condition (1,2), polarity (plus,minus)}
                        SpikeTrains={SpikeTrainsA_plus,SpikeTrainsA_minus;SpikeTrainsB_plus,SpikeTrainsB_minus};

                        % specify params to be used
                        clear paramsIN
                        paramsIN.durA_msec=dur_sec*1000;
                        paramsIN.durB_msec=dur_sec*1000;
                        paramsIN.CF_Hz=CF_kHz(Fiber_Number)*1000;

                        [SACSCCfunctions,SACSCCmetrics,paramsOUT] = CCCanal(SpikeTrains,paramsIN,0);

                        if i==1 % aided+impaired
                            Difcor(phone_index,level_index,gain_index,Fiber_Number,spont_index) = SACSCCmetrics.DCpeak_A;
                            Sumcor(phone_index,level_index,gain_index,Fiber_Number,spont_index) = SACSCCmetrics.SCpeak_A;
                            Env(phone_index,level_index,gain_index,Fiber_Number,spont_index) = SACSCCmetrics.CCCenv;
                            Tfs(phone_index,level_index,gain_index,Fiber_Number,spont_index) = SACSCCmetrics.CCCtfs;
                        elseif i==2 % normal
                            Difcor_Normal(phone_index,level_index,Fiber_Number,spont_index) = SACSCCmetrics.DCpeak_A;
                            Sumcor_Normal(phone_index,level_index,Fiber_Number,spont_index) = SACSCCmetrics.SCpeak_A;
                            Env_Normal(phone_index,level_index,Fiber_Number,spont_index) = SACSCCmetrics.CCCenv;
                            Tfs_Normal(phone_index,level_index,Fiber_Number,spont_index) = SACSCCmetrics.CCCtfs;
                        end
                        save NAG PathName FileName Rate Rate_Normal Difcor Difcor_Normal Sumcor Sumcor_Normal Env Env_Normal Tfs Tfs_Normal phones levels gains sponts;
                    end % end numreps (impaired & normal)

                    spont_index=spont_index+1;    
                end % end SR_sps
            end % end Fiber_Number
            gain_index=gain_index+1;
        end % end Gain_Adjust
        level_index=level_index+1;
    end % end OALevel_dBSPL
    phone_index = phone_index+1;
end % end phone


