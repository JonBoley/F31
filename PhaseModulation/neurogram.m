function Neurogram = neurogram(input,Fs,OALevel_dBSPL,lowCF,highCF,numCFs,impairment,aided)

% Check for Neural Model
if strcmp(which('catmodel_IHC'),'') % if we cannot find this function
    pathtool; % Let user add the appropriate directory to the path
    errordlg('Cannot Find Neural Model!');
    return;
end

if nargin<2
    % Get Input Audio File Name
    [FileName,PathName,FilterIndex] = uigetfile('*.wav',...
        'Select Input Audio File','heed_100.wav');
    pause(0.1); % just wait for window to disappear
    [input,Fs,NBITS]=wavread([PathName FileName]);
end
if nargin<3
    OALevel_dBSPL = 70;
end
if nargin <4
    lowCF=80; %Hz
end
if nargin <5
    highCF=4000; %Hz
end
if nargin <6
    numCFs=60;
end

% Neural Model Params
spread=log2(highCF/lowCF); %number of octaves to span
CF_kHz=(lowCF*2.^(0:spread/(numCFs-1):spread))/1000;
Cohc=1.0*ones(length(CF_kHz),1);  % initial OHC health
Cihc=1.0*ones(length(CF_kHz),1);  % initial IHC health
Nreps=10; % number of repetitions
ANmodel_Fs_Hz=100e3; % 100kHz sampling

% Model hearing loss
FREQUENCIES = [250,500,1000,2000,3000,4000,6000]; % Hz
dBLoss = [20 20 30 40 45 50 50]; % Mild Hearing Loss (according to Bruce (ISAAR 2007)
if nargin==7
%     impairment = 'mixed';
    switch impairment
        case 'mixed'
            Dsd_OHC_Loss = 2/3*dBLoss;%zeros(1,length(FREQUENCIES));%24 *ones(1,length(FREQUENCIES));
        case 'ohc'
            Dsd_OHC_Loss = dBLoss;
        otherwise %case ihc
            Dsd_OHC_Loss = 0*dBLoss;
    end
    [Cohc,Cihc,OHC_Loss]=fitaudiogram(FREQUENCIES,dBLoss,Dsd_OHC_Loss);
    % Interpolate impairment
    Cihc = interp1(FREQUENCIES,Cihc,CF_kHz*1000,'linear',min(Cihc));
    Cohc = interp1(FREQUENCIES,Cohc,CF_kHz*1000,'linear',min(Cohc));
end
if ((nargin==8) && (aided==true))
    % NAL equations
    % frequency shaping (250,500,1000,2000,3000,4000,6000Hz)
    k_NAL = [-18 -8 1 -1 -2 -2 -2]; % dB 
    H_3FA = sum(dBLoss(2:4)) / 3; % sum up loss at 500Hz,1kHz,2kHz
    X = 0.15*H_3FA;
    R = 0.31; % NAL-RP formula  (not quite half-gain rule)
    NAL_IG = X + R.*dBLoss + k_NAL; % insertion gain
    NAL_IG = max(0,NAL_IG); % no negative gain
    
    NAL = max(0,interp1(FREQUENCIES,NAL_IG,CF_kHz*1000));
else
    NAL = zeros(length(CF_kHz),1);
end


%% Analyze the data
input = input/max(input); % normalize
dur_sec = length(input)/Fs;

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
    error('RESAMPLING CHANGED input dBSPL by %f dB',dBSPL_after-dBSPL_before)
end

SynOut = cell(numCFs,3); %(3 populations of spont rates)
w_spont = [.6 .2 .2]; %[60%HSR 20%MSR 20%LSR]

for Fiber_Number=1:numCFs
    fprintf('.'); if (mod(Fiber_Number,5)==0), fprintf(' '); end
    if (Fiber_Number==numCFs), fprintf('\n'); end
    
    adjustment = max(OALevel_dBSPL+NAL(Fiber_Number)-dBSPL_after,-dBSPL_after); % dB
    input_model_fiber = input_model*10^((adjustment)/20);
    input_model_fiber = refit_waveform(input_model_fiber,ANmodel_Fs_Hz,dur_sec*1000);
    input_model_fiber = window_waveform(input_model_fiber,ANmodel_Fs_Hz,dur_sec*1000);

    spont_index=1;
    for fibertype = 3:-1:1 %High, Med, Low Spont

        % vihc = catmodel_IHC(pin,CF,nrep,binwidth,reptime,cohc,cihc);
        % [synout,psth] = catmodel_Synapse(vihc,CF,nrep,binwidth,reptime,fibertype,implnt);
        vihc = catmodel_IHC(input_model_fiber.',CF_kHz(Fiber_Number)*1000,Nreps,1/ANmodel_Fs_Hz,dur_sec+0.200,Cohc(Fiber_Number),Cihc(Fiber_Number));
        [SynOut{Fiber_Number,spont_index},psth] = catmodel_Synapse(vihc,CF_kHz(Fiber_Number)*1000,Nreps,1/ANmodel_Fs_Hz,dur_sec+0.200,fibertype,0);

        if Fiber_Number==1 && spont_index==1
            Neurogram = zeros(length(SynOut{1,1}),numCFs);
        end
        if ~isempty(SynOut{Fiber_Number,spont_index})
            Neurogram(:,Fiber_Number) = Neurogram(:,Fiber_Number) + ...
                w_spont(spont_index)*SynOut{Fiber_Number,spont_index}';
        end
        spont_index=spont_index+1;
    end % end SR_sps
    
%     N_win = round(.008*ANmodel_Fs_Hz); % 8ms
%     W = hann(N_win)/sum(hann(N_win_short));
%     Neurogram(:,Fiber_Number) = filter(W,1,Neurogram(:,Fiber_Number));
    
end % end Fiber_Number

Neurogram = Neurogram';
if nargout<1
    scrsz = get(0,'ScreenSize');
    figure('Position',[scrsz(3)/2-300 scrsz(4)/2-100 600 200]);
    imagesc(Neurogram); axis off;
end

end
