% function NAG(analysis)
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

EMAIL_NOTIFICATION = 1; %turn on/off email notification

GEN_SPIKES = 1;
ANALYZE_SPIKES = 1;

strategy_list = {'short','avg','rate','env','tfs'};
STRATEGY = 4;

START_PHONE = 4;

levels = 65;%[45 65 85];
gains = -40:5:40;
note = 'Sep_07_11';%datestr(now,'mmm_dd_yy'); %attach note to end of file name

NumLabs = 2; %use zero for max number of parallel processors
if NumLabs
    if (NumLabs>1 && NumLabs~=matlabpool('size'))
        if matlabpool('size'), matlabpool close; end %close open labs
        eval(sprintf('matlabpool open %d;',NumLabs));
    end
else
    if matlabpool('size')<8 %if we're not already maxed out
        if matlabpool('size'), matlabpool close; end %close open labs
        matlabpool; % enable parallel computing with max cores
    end
end

%% Get username/password for email notification
if(EMAIL_NOTIFICATION)
    [login password] = logindlg('Title','GMail Login');
end

%% Initialization Routine
% Get Input Audio File Name
PathName = 'TIMIT\test\dr1\faks0\';
FileName = 'sa1.wav';
if (true) % change to false if you want to ask the user for an input file
    phonemeindx = textread([PathName FileName(1:end-3) 'phn'],'%*d %d %*s');
else
    [FileName,PathName,FilterIndex] = uigetfile('*.wav',...
        'Select Input Audio File',[PathName 'sa1.wav']);
    % Get TIMIT phoneme structure if needed
    reply = menu('Load TIMIT phoneme structure?','YES','NO');
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

% Neural Model Params
lowCF=250; %Hz
highCF=4000; %Hz
numCFs=30;
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
impairment = 'mixed';
switch impairment
    case 'mixed'
        Dsd_OHC_Loss = 2/3*dBLoss;%zeros(1,length(FREQUENCIES));%24 *ones(1,length(FREQUENCIES));
    case 'ohc'
        Dsd_OHC_Loss = dBLoss;
    otherwise %case ihc
        Dsd_OHC_Loss = 0*dBLoss;
end
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
fprintf('Neural Analysis of Gain\nby Jon Boley\n');
fprintf('jdboley(a)purdue.edu\n%s\n',file.date);


%% Analyze the data
[input_orig,Fs] = readnist([PathName FileName]);
input_orig = input_orig/max(input_orig);

% Apply NAL Gain
NAL_filter_freqs = 2.^interp1(0.5:6.5,log2(FREQUENCIES),0:7,'linear','extrap')/(Fs/2);
NAL_filter_gains = 10.^(interp1(0.5:6.5,NAL_IG,0:7,'linear','extrap')/20);
b = firpm(16,NAL_filter_freqs,NAL_filter_gains);

phones = START_PHONE:length(phonemeindx);

sponts = [50, 5, 0.25];
N_win_short = round(.000256*ANmodel_Fs_Hz); % 256us
N_win_long = round(.008*ANmodel_Fs_Hz); % 8ms
W_win_short = hann(N_win_short)/sum(hann(N_win_short));
W_win_long = hann(N_win_long)/sum(hann(N_win_long));
w_spont = [.6 .2 .2]; %[60%HSR 20%MSR 20%LSR]

level_index=1;
for OALevel_dBSPL=levels
    disp(sprintf('Processing level of %ddBSPL',OALevel_dBSPL));
    
    PhonemeLevel = NaN*ones(length(phones),length(levels));
    OptimumGain = NaN*ones(length(phones),length(levels));

    dBSPL_before=20*log10(sqrt(mean(input_orig.^2))/(20e-6));
    input_orig = input_orig*10^((OALevel_dBSPL-dBSPL_before)/20);
    
    % Apply NAL gain
    input_orig_NAL = filter(b,1,input_orig);
    % Remove group delay (to maintain phoneme alignment)
    % (NOTE: it would be interesting to play with the group delay here)
    input_orig_NAL = [input_orig_NAL(ceil(length(b)/2):end); zeros(floor(length(b)/2),1)];
    
    phone_index=min(phones);
    for phone=phones
        disp(sprintf('Processing phoneme #%d',phone));

        %%%%%%% Will need to include all previous phones
        % Note: neurogram for each phoneme has been extended in time
        input=input_orig(1:min(phonemeindx(phone),length(input_orig)));
        input_NAL=input_orig_NAL(1:min(phonemeindx(phone),length(input_orig_NAL)));
        dur_sec = length(input)/Fs;
        
        
        if phone_index<=1
            PhonemeLevel(phone_index,level_index) = ...
                20*log10(sqrt(mean(input(1:min(phonemeindx(phone),length(input))).^2))/(20e-6));
        else % if phone_index>1
            PhonemeLevel(phone_index,level_index) = ...
                20*log10(sqrt(mean(input(phonemeindx(phone-1):min(phonemeindx(phone),length(input))).^2))/(20e-6));
            
            for phone_index2=2:phone
                % Determine OptimumGain for each phone
                [short,avg,rate,env,tfs] = OptimalGain2('archive\',levels,impairment,strategy_list{STRATEGY},phone_index2-1,gains,note);
                eval(sprintf('OptimumGain(phone_index2-1,:)=%s;',strategy_list{STRATEGY}));
                fprintf('Using optimal gain (%s) of %ddB on phone %d\n',strategy_list{STRATEGY},OptimumGain(phone_index2-1,:),phone_index2-1);
            end
            
            % Apply OptimumGain to previous phones. (NOTE: We could add attack/release here)
            StartIndex = 1; % start of phone
            EndIndex = min(phonemeindx(1),length(input_NAL)); % end of phone
            input_NAL(StartIndex:EndIndex) = ...
                input_NAL(StartIndex:EndIndex)*10^((OptimumGain(1,level_index))/20);
            for PhoneNum=2:(phone-1)
                StartIndex = phonemeindx(PhoneNum-1)+1;
                EndIndex = min(phonemeindx(PhoneNum),length(input_NAL));
                input_NAL(StartIndex:EndIndex) = ...
                    input_NAL(StartIndex:EndIndex)*10^((OptimumGain(PhoneNum,level_index))/20);
            end
        end

        % Resample original version for AN model
        dBSPL_before=20*log10(sqrt(mean(input.^2))/(20e-6));
        sfreq=Fs;
        sfreqNEW=ANmodel_Fs_Hz;
        P=round(sfreqNEW/10); Q=round(sfreq/10);  %Integers used to up sample
        if(P/Q*sfreq~=sfreqNEW), disp('Integer sfreq conversion NOT exact'); end
        Nfir=30;  % proportional to FIR filter length used for resampling: higher Nfir, better accuracy & longer comp time
        input_model_normal=resample(input,P,Q,Nfir);
        dBSPL_after=20*log10(sqrt(mean(input_model_normal.^2))/(20e-6));
        if abs(dBSPL_before-dBSPL_after)>2;
            error('RESAMPLING CHANGED input dBSPL by %f dB',dBSPL_after-dBSPL_before)
        end
        input_model_normal = refit_waveform(input_model_normal,ANmodel_Fs_Hz,dur_sec*1000);
        input_model_normal = window_waveform(input_model_normal,ANmodel_Fs_Hz,dur_sec*1000);
        
        % Resample NAL version for AN model
        dBSPL_before=20*log10(sqrt(mean(input_NAL.^2))/(20e-6));
        sfreq=Fs;
        sfreqNEW=ANmodel_Fs_Hz;
        P=round(sfreqNEW/10); Q=round(sfreq/10);  %Integers used to up sample
        if(P/Q*sfreq~=sfreqNEW), disp('Integer sfreq conversion NOT exact'); end
        Nfir=30;  % proportional to FIR filter length used for resampling: higher Nfir, better accuracy & longer comp time
        input_model_NAL=resample(input_NAL,P,Q,Nfir);
        dBSPL_after=20*log10(sqrt(mean(input_model_NAL.^2))/(20e-6));
        if abs(dBSPL_before-dBSPL_after)>2;
            error('RESAMPLING CHANGED input dBSPL by %f dB',dBSPL_after-dBSPL_before)
        end
        input_model_NAL = refit_waveform(input_model_NAL,ANmodel_Fs_Hz,dur_sec*1000);
        input_model_NAL = window_waveform(input_model_NAL,ANmodel_Fs_Hz,dur_sec*1000);
        
        if phone_index<=1
            StartIndex_mdl = 1;
        else
            StartIndex_mdl = round(phonemeindx(phone_index-1)/Fs*ANmodel_Fs_Hz+1);
        end
        EndIndex_mdl = min(round(phonemeindx(phone_index)/Fs*ANmodel_Fs_Hz),length(input_model_NAL));
        
        Mydiff1 = zeros(length(levels),length(gains));
        Mydiff2 = zeros(length(levels),length(gains));
        
        gain_index=1;
        for Gain_Adjust=gains
            disp(sprintf('Adjusting prescribed gain by %ddB (%s)',Gain_Adjust,datestr(now)));
            
            Gain_Adjust_array = [ones(StartIndex_mdl-1,1); 10^(Gain_Adjust/20)*ones(EndIndex_mdl-StartIndex_mdl+1,1)];
            input_model_NAL2 = input_model_NAL .* Gain_Adjust_array;
            input_model_normal2 = input_model_normal .* Gain_Adjust_array;

            % initialization
            SynOutA = cell(numCFs,1);
            for z=1:numCFs, SynOutA{z} = cell(length(sponts),1); end
            SynOutB = cell(numCFs,1);
            for z=1:numCFs, SynOutB{z} = cell(length(sponts),1); end
            if (GEN_SPIKES)
                SpikesA_plus = cell(numCFs,1);
                for z=1:numCFs, SpikesA_plus{z} = cell(length(sponts),1); end
                SpikesA_minus = cell(numCFs,1);
                for z=1:numCFs, SpikesA_minus{z} = cell(length(sponts),1); end
                SpikesB_plus = cell(numCFs,1);
                for z=1:numCFs, SpikesB_plus{z} = cell(length(sponts),1); end
                SpikesB_minus = cell(numCFs,1);
                for z=1:numCFs, SpikesB_minus{z} = cell(length(sponts),1); end
            end
            if (ANALYZE_SPIKES)
                Rate   = cell(numCFs,1);
                for z=1:numCFs, Rate{z} = NaN*ones(length(sponts),1); end
                Difcor = cell(numCFs,1);
                for z=1:numCFs, Difcor{z} = NaN*ones(length(sponts),1); end
                Sumcor = cell(numCFs,1);
                for z=1:numCFs, Sumcor{z} = NaN*ones(length(sponts),1); end
                Env = cell(numCFs,1);
                for z=1:numCFs, Env{z} = NaN*ones(length(sponts),1); end
                Tfs = cell(numCFs,1);
                for z=1:numCFs, Tfs{z} = NaN*ones(length(sponts),1); end
            end

            %default PARAMS
            if (ANALYZE_SPIKES)
                paramsIN = cell(numCFs,1);
            end
            
            len_Neurogram = EndIndex_mdl-StartIndex_mdl+1;
            neurogramA1 = zeros(len_Neurogram,numCFs);%impaired
            neurogramA2 = zeros(len_Neurogram,numCFs);%impaired
            neurogramA3 = zeros(len_Neurogram,numCFs);%impaired
            neurogramB1 = zeros(len_Neurogram,numCFs);%normal
            neurogramB2 = zeros(len_Neurogram,numCFs);%normal
            neurogramB3 = zeros(len_Neurogram,numCFs);%normal
            
            for Fiber_Number=1:numCFs
                fprintf('.'); if (mod(Fiber_Number,5)==0), fprintf(' '); end
                if (Fiber_Number==numCFs), fprintf('\n'); end
                %disp(sprintf('Processing fiber #%d of %d',Fiber_Number,numCFs));
                
                if (ANALYZE_SPIKES)
                    % specify params to be used
                    paramsIN{Fiber_Number}.durA_msec=dur_sec*1000;
                    paramsIN{Fiber_Number}.durB_msec=dur_sec*1000;
                    paramsIN{Fiber_Number}.CF_A_Hz=CF_kHz(Fiber_Number)*1000;
                    paramsIN{Fiber_Number}.CF_B_Hz=CF_kHz(Fiber_Number)*1000;
                    paramsIN{Fiber_Number}.MAXspikes=3000;
                    paramsIN{Fiber_Number}.PSD_LHfreqs_Hz=[0 64; 0 50];  %additional freq ranges to compute CCCenv for
                end
                
                spont_index=1;
                for fibertype = 3:-1:1 %High, Med, Low Spont
                    %                     SR_sps = sponts
                    %disp(sprintf('Processing spont rate of %d',SR_sps));
                    % Run fiber (A+)
                    %                     [timeout,meout,c1filterout,c2filterout,c1vihc,c2vihc,vihc,sout,psth500k] ...
                    %                         = zbcatmodel(input_model_NAL.',CF_kHz(Fiber_Number)*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.05,Cohc_impaired,Cihc_impaired,SR_sps);
                    vihc = catmodel_IHC(input_model_NAL2.',CF_kHz(Fiber_Number)*1000,1,1/ANmodel_Fs_Hz,dur_sec+1.00,Cohc_impaired,Cihc_impaired);
                    [sout,psth] = catmodel_Synapse(vihc,CF_kHz(Fiber_Number)*1000,1,1/ANmodel_Fs_Hz,dur_sec+1.00,fibertype,1);
                    if (GEN_SPIKES)
                        [sptimes nspikes]= SGfast([1/ANmodel_Fs_Hz, Nreps],sout);
                        sptimes=sptimes(sptimes>=StartIndex_mdl/ANmodel_Fs_Hz & sptimes<=EndIndex_mdl/ANmodel_Fs_Hz); 
                        nspikes=length(sptimes);
                        Nreps2=Nreps; 
                        while nspikes<2000, 
                            Nreps2=Nreps2*2; 
                            if (Nreps2>2000 && Nreps2>=1), sptimes=0; break; end; 
                            [sptimes nspikes]= SGfast([1/ANmodel_Fs_Hz, Nreps2],sout); 
                            sptimes=sptimes(sptimes>=StartIndex_mdl/ANmodel_Fs_Hz & sptimes<=EndIndex_mdl/ANmodel_Fs_Hz);
                            nspikes=length(sptimes);
                        end
                        % Format spikes into NEL spikes format then cell array
                        NELspikes=ANmodelSTs2nel(sptimes,Nreps2);
                        SpikeTrainsA_plus=nelSTs2cell(NELspikes);
                        if length(SpikeTrainsA_plus)<Nreps2
                            SpikeTrainsA_plus = [SpikeTrainsA_plus cell(1,Nreps2-length(SpikeTrainsA_plus))];
                        end
                        SpikesA_plus{Fiber_Number}{spont_index} = SpikeTrainsA_plus;
                    end
                    sout=sout(StartIndex_mdl:EndIndex_mdl);
                    SynOutA{Fiber_Number}{spont_index} = sout;
                    len_SynOutA = length(sout);
                    if (ANALYZE_SPIKES)
                        Rate{Fiber_Number}(spont_index) = mean(sout); % SAVE RATE INFO
                    end

                    % Run fiber (A-)
                    %                     [timeout,meout,c1filterout,c2filterout,c1vihc,c2vihc,vihc,sout,psth500k] ...
                    %                         = zbcatmodel(-input_model_NAL.',CF_kHz(Fiber_Number)*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.05,Cohc_impaired,Cihc_impaired,SR_sps);
                    vihc = catmodel_IHC(-input_model_NAL2.',CF_kHz(Fiber_Number)*1000,1,1/ANmodel_Fs_Hz,dur_sec+1.00,Cohc_impaired,Cihc_impaired);
                    [sout,psth] = catmodel_Synapse(vihc,CF_kHz(Fiber_Number)*1000,1,1/ANmodel_Fs_Hz,dur_sec+1.00,fibertype,1);
                    if (GEN_SPIKES)
                        [sptimes nspikes]= SGfast([1/ANmodel_Fs_Hz, Nreps],sout);
                        sptimes=sptimes(sptimes>=StartIndex_mdl/ANmodel_Fs_Hz & sptimes<=EndIndex_mdl/ANmodel_Fs_Hz); 
                        nspikes=length(sptimes);
                        Nreps2=Nreps; 
                        while nspikes<2000, 
                            Nreps2=Nreps2*2; 
                            if (Nreps2>2000 && Nreps2>=1), sptimes=0; break; end; 
                            [sptimes nspikes]= SGfast([1/ANmodel_Fs_Hz, Nreps2],sout); 
                            sptimes=sptimes(sptimes>=StartIndex_mdl/ANmodel_Fs_Hz & sptimes<=EndIndex_mdl/ANmodel_Fs_Hz);
                            nspikes=length(sptimes);
                        end
                        % Format spikes into NEL spikes format then cell array
                        NELspikes=ANmodelSTs2nel(sptimes,Nreps2);
                        SpikeTrainsA_minus=nelSTs2cell(NELspikes);
                        if length(SpikeTrainsA_minus)<Nreps2
                            SpikeTrainsA_minus = [SpikeTrainsA_minus cell(1,Nreps2-length(SpikeTrainsA_minus))];
                        end
                        SpikesA_minus{Fiber_Number}{spont_index} = SpikeTrainsA_minus;
                    end

                    % Run fiber (B+) - always normal
                    %                     [timeout,meout,c1filterout,c2filterout,c1vihc,c2vihc,vihc,sout,psth500k] ...
                    %                         = zbcatmodel(input_model_normal.',CF_kHz(Fiber_Number)*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.05,Cohc,Cihc,SR_sps);
                    vihc = catmodel_IHC(input_model_normal2.',CF_kHz(Fiber_Number)*1000,1,1/ANmodel_Fs_Hz,dur_sec+1.00,Cohc,Cihc);
                    [sout,psth] = catmodel_Synapse(vihc,CF_kHz(Fiber_Number)*1000,1,1/ANmodel_Fs_Hz,dur_sec+1.00,fibertype,1);
                    if (GEN_SPIKES)
                        [sptimes nspikes]= SGfast([1/ANmodel_Fs_Hz, Nreps],sout);
                        sptimes=sptimes(sptimes>=StartIndex_mdl/ANmodel_Fs_Hz & sptimes<=EndIndex_mdl/ANmodel_Fs_Hz); 
                        nspikes=length(sptimes);
                        Nreps2=Nreps; 
                        while nspikes<2000, 
                            Nreps2=Nreps2*2; 
                            if (Nreps2>2000 && Nreps2>=1), sptimes=0; break; end; 
                            [sptimes nspikes]= SGfast([1/ANmodel_Fs_Hz, Nreps2],sout); 
                            sptimes=sptimes(sptimes>=StartIndex_mdl/ANmodel_Fs_Hz & sptimes<=EndIndex_mdl/ANmodel_Fs_Hz);
                            nspikes=length(sptimes);
                        end
                        % Format spikes into NEL spikes format then cell array
                        NELspikes=ANmodelSTs2nel(sptimes,Nreps2);
                        SpikeTrainsB_plus=nelSTs2cell(NELspikes);
                        if length(SpikeTrainsB_plus)<Nreps2
                            SpikeTrainsB_plus = [SpikeTrainsB_plus cell(1,Nreps2-length(SpikeTrainsB_plus))];
                        end
                        SpikesB_plus{Fiber_Number}{spont_index} = SpikeTrainsB_plus;
                    end
                    sout=sout(StartIndex_mdl:EndIndex_mdl);
                    SynOutB{Fiber_Number}{spont_index} = sout;
                    len_SynOutB = length(sout);
                    
                    % Run fiber (B-) - always normal
                    %                     [timeout,meout,c1filterout,c2filterout,c1vihc,c2vihc,vihc,sout,psth500k] ...
                    %                         = zbcatmodel(-input_model_normal.',CF_kHz(Fiber_Number)*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.05,Cohc,Cihc,SR_sps);
                    vihc = catmodel_IHC(-input_model_normal2.',CF_kHz(Fiber_Number)*1000,1,1/ANmodel_Fs_Hz,dur_sec+1.00,Cohc,Cihc);
                    [sout,psth] = catmodel_Synapse(vihc,CF_kHz(Fiber_Number)*1000,1,1/ANmodel_Fs_Hz,dur_sec+1.00,fibertype,1);
                    if (GEN_SPIKES)
                        [sptimes nspikes]= SGfast([1/ANmodel_Fs_Hz, Nreps],sout);
                        sptimes=sptimes(sptimes>=StartIndex_mdl/ANmodel_Fs_Hz & sptimes<=EndIndex_mdl/ANmodel_Fs_Hz); 
                        nspikes=length(sptimes);
                        Nreps2=Nreps; 
                        while nspikes<2000, 
                            Nreps2=Nreps2*2; 
                            if (Nreps2>2000 && Nreps2>=1), sptimes=0; break; end; 
                            [sptimes nspikes]= SGfast([1/ANmodel_Fs_Hz, Nreps2],sout); 
                            sptimes=sptimes(sptimes>=StartIndex_mdl/ANmodel_Fs_Hz & sptimes<=EndIndex_mdl/ANmodel_Fs_Hz);
                            nspikes=length(sptimes);
                        end
                        % Format spikes into NEL spikes format then cell array
                        NELspikes=ANmodelSTs2nel(sptimes,Nreps2);
                        SpikeTrainsB_minus=nelSTs2cell(NELspikes);
                        if length(SpikeTrainsB_minus)<Nreps2
                            SpikeTrainsB_minus = [SpikeTrainsB_minus cell(1,Nreps2-length(SpikeTrainsB_minus))];
                        end
                        SpikesB_minus{Fiber_Number}{spont_index} = SpikeTrainsB_minus;
                    end
                    if (ANALYZE_SPIKES)
                        % Calculate env & tfs coding
                        % Spike Analysis
                        % Organize variables for CCCanal
                        SpikeTrains={SpikeTrainsA_plus,SpikeTrainsA_minus;SpikeTrainsB_plus,SpikeTrainsB_minus};

                        try
                            %	[SACSCCfunctions,SACSCCmetrics,paramsOUT] = CCCanal(SpikeTrains,paramsIN,0);
                            [SACSCCfunctions,SACSCCmetrics,paramsOUT] = CCCanal_4(SpikeTrains,paramsIN{Fiber_Number},0);
                        catch exception
                            dbstop;
                            rethrow(exception);
                        end

                        Difcor{Fiber_Number}(spont_index) = SACSCCmetrics.DCpeak_A;
                        Sumcor{Fiber_Number}(spont_index) = ...
                            SACSCCmetrics.SCpeaks_A(find(strcmp('IFFTraw',SACSCCmetrics.SCpeaks_legend)));
                        Env{Fiber_Number}(spont_index) = ...
                            SACSCCmetrics.CCCenvs(find(strcmp('0-300, subBIAS',SACSCCmetrics.CCCenvs_legend)));
                        Tfs{Fiber_Number}(spont_index) = SACSCCmetrics.CCCtfs;
                    end %if

                    if ~isempty(SynOutA{Fiber_Number}{spont_index})
                        neurogramA1(:,Fiber_Number) = neurogramA1(:,Fiber_Number) + ...
                            w_spont(spont_index)*SynOutA{Fiber_Number}{spont_index}';
                        neurogramB1(:,Fiber_Number) = neurogramB1(:,Fiber_Number) + ...
                            w_spont(spont_index)*SynOutB{Fiber_Number}{spont_index}';
                    end
                    spont_index=spont_index+1;
                end % end fibertype (or SR_sps with the old model)
                neurogramA3(:,Fiber_Number) = neurogramA1(:,Fiber_Number);
                neurogramB3(:,Fiber_Number) = neurogramB1(:,Fiber_Number); % unfiltered version
                neurogramA2(:,Fiber_Number) = filter(W_win_long,1,neurogramA1(:,Fiber_Number));
                neurogramB2(:,Fiber_Number) = filter(W_win_long,1,neurogramB1(:,Fiber_Number)); % apply long window
                neurogramA1(:,Fiber_Number) = filter(W_win_short,1,neurogramA1(:,Fiber_Number));
                neurogramB1(:,Fiber_Number) = filter(W_win_short,1,neurogramB1(:,Fiber_Number)); % apply short window
                
            end % end Fiber_Number
%             Mydiff1(level_index,gain_index) = mean(mean(abs(neurogramB1-neurogramA1)));
%             Mydiff2(level_index,gain_index) = mean(mean(abs(neurogramB2-neurogramA2)));

            % Save file
            MatFileName = sprintf('archive\\%1.0fdBSPL\\%s\\%s\\phone%d\\%1.0fdBgain_%s',...
                OALevel_dBSPL,impairment,strategy_list{STRATEGY},phone,Gain_Adjust,note);
            if ~exist(sprintf('archive\\%1.0fdBSPL\\%s\\%s\\phone%d',OALevel_dBSPL,impairment,strategy_list{STRATEGY},phone),'dir')
                mkdir(sprintf('archive\\%1.0fdBSPL\\%s\\%s\\phone%d',OALevel_dBSPL,impairment,strategy_list{STRATEGY},phone));
            end

            save(MatFileName,'PathName','FileName','neurogramA1','neurogramA2','neurogramA3');
            save(MatFileName,'PhonemeLevel','OptimumGain','-append');
            if (Gain_Adjust==0)
                save(MatFileName,'neurogramB1','neurogramB2','neurogramB3','-append');
            end
            if (GEN_SPIKES)
                save(MatFileName,'SpikesA_plus','SpikesA_minus','-append');
                if Gain_Adjust==0,
                    save(MatFileName,'SpikesB_plus','SpikesB_minus','-append');
                end
            end
            if (ANALYZE_SPIKES)
                save(MatFileName,'Difcor','Sumcor','Env','Tfs','-append');
            end

            clear SynOutA SynOutB
            if (GEN_SPIKES)
                clear SpikesA_plus SpikesA_minus SpikesB_plus SpikesB_minus
            end
            if (ANALYZE_SPIKES)
                clear Rate Difcor Sumcor Env Tfs
            end

            if(EMAIL_NOTIFICATION)
                try
                    email_text = ...
                        sprintf('Calculation completed (%s):\nLevel: %d\nPhone: %d\nGain: %d\n\nSent %s\n',...
                        strategy_list{STRATEGY},OALevel_dBSPL,phone,Gain_Adjust,datestr(now));
                    EmailNotification(login,password,'jdboley@purdue.edu','Matlab Update',email_text);
                    fprintf('Notification Email Sent\n');
                catch
                    fprintf('[Error Sending Notification Email]\n');
                end
            end
            
            gain_index=gain_index+1;
        end % end Gain_Adjust
        
        phone_index = phone_index+1;
    end % end phone

	level_index=level_index+1;
end % end OALevel_dBSPL


