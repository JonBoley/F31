function STMPsave_STs_1unitNstim_interleaved(ExpDate,UnitName,STIMtype)
% 6/5/09: Added functionality from BBNsave_STs_1unitNstim( )
% FROM: STMPsave_STs_1unitNstim_EHrBFi(ExpDate,UnitName)
% FROM: STMPsave_STs_1unitNstim_EHvNrBFi(ExpDate,UnitName)
% M. Heinz Jun 11, 2007
% Processes 1 unit for STMP analysis, for which data collected was from 1
% unit and multiple stimuli (typically varied sampling rates)
%
% Stimuli: EHrBFi, EHvNBFi, TrBFi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global STMP_dir STMP_ExpList
global FeaturesText FormsAtHarmonicsText InvertPolarityText

STIMtypes_allowed={'EHrBFi','EHvNrBFi','TrBFi','BBN_SX','BBNrBFi','WAVreBFi'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Verify parameters and experiment, unit are valid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TESTunitNUM=4;  % 1: ARO (111804, 1.28), 2: ISH (041805, 2.04), 3: Purdue1 (111606, 2.09); 4: Purdue2 (041306, 1.03),
%%%% Specify ExpDate if not provided
if ~exist('ExpDate','var')
    if TESTunitNUM==1
        ExpDate='111804'; UnitName='1.28'; STIMtype='EHrBFi';
    elseif TESTunitNUM==2
        ExpDate='041805'; UnitName='2.04'; STIMtype='EHvNrBFi';
    elseif TESTunitNUM==3
        ExpDate='111606'; UnitName='2.09'; STIMtype='EHvNrBFi';  % Purdue #1
    elseif TESTunitNUM==4
        ExpDate='041307'; UnitName='2.03'; STIMtype='EHvNrBFi';  % Purdue #2
    elseif TESTunitNUM==5
        ExpDate='041805'; UnitName='2.04'; STIMtype='TrBFi';   % TONE re BFI, with Tone RLV
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CHECK STIMtype is setup!
if ~sum(strcmp(STIMtypes_allowed,STIMtype))
    beep
    error(sprintf('STIMtype: %s NOT SETUP YET!!',STIMtype))
end

%%%% Find the full Experiment Name
ExpDateText=strcat('20',ExpDate(end-1:end),'_',ExpDate(1:2),'_',ExpDate(3:4));
for i=1:length(STMP_ExpList)
    if ~isempty(strfind(STMP_ExpList{i},ExpDateText))
        ExpName=STMP_ExpList{i};
        break;
    end
end
if ~exist('ExpName','var')
    disp(sprintf('***ERROR***:  Experiment: %s not found\n   Experiment List:',ExpDate))
    disp(strvcat(STMP_ExpList))
    beep
    error('STOPPED');
end

%%%% Parse out the Track and Unit Number
TrackNum=str2num(UnitName(1:strfind(UnitName,'.')-1));
UnitNum=str2num(UnitName(strfind(UnitName,'.')+1:end));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Setup related directories
RAWdata_dir=fullfile(STMP_dir,'ExpData',ExpName);
eval(['cd ''' RAWdata_dir ''''])
UNITinfo_dir=fullfile(STMP_dir,'ExpData',ExpName,'UNITinfo');   % For general unit info (BF, SR, bad lines, ...)
if ~exist('UNITinfo','dir')
    mkdir('UNITinfo');
end
STMPanal_dir=fullfile(STMP_dir,'ExpData',ExpName,'STMPanalyses');   % For STMP analyses (Spike Trains, PSTs, PerHist, DFTs, SAC/SCCs, ...)
if ~exist('STMPanalyses','dir')
    mkdir('STMPanalyses');
end
disp(sprintf('... STMP: Saving ORIGINAL "%s" SpikeTrains for:  Experiment: ''%s''; Unit: %d.%02d',STIMtype,ExpName,TrackNum,UnitNum))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Verify unitINFO exists, if not create
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
unitINFO_filename=sprintf('unitINFO.%d.%02d.mat',TrackNum,UnitNum);
eval(['ddd=dir(''' fullfile(UNITinfo_dir,unitINFO_filename) ''');'])
if isempty(ddd)
    eval(['cd ''' RAWdata_dir ''''])
    STMPsave_UnitINFO(TrackNum,UnitNum,1);  % skip verify
end
eval(['load ''' fullfile(UNITinfo_dir,unitINFO_filename) ''''])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% See if STs_1unitNstim exists already
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
STs_1unitNstim_filename=sprintf('STs_1unitNstim.%d.%02d.%s.mat',TrackNum,UnitNum,STIMtype);
eval(['ddd=dir(''' fullfile(STMPanal_dir,STs_1unitNstim_filename) ''');'])
if ~isempty(ddd)
    beep
    TEMP = input(sprintf('File: ''%s'' already exists!!\n  ***** Do you want to re-run STs_1unitNstim, or leave it as is?  [0]: LEAVE AS IS; RERUN: 1;  ',STs_1unitNstim_filename));
    if isempty(TEMP)
        TEMP=0;
    end
    if TEMP~=1
        beep
        disp(sprintf(' FILE NOT ALTERED\n'));
        return;
    else
        disp(' ... Re-Running STs_1unitNstim - SAVING NEW FILE!');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup STs_1unitNstim data structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save basic INFO (see unitINFO for rest)
STs_1unitNstim.Info.ExpName=unitINFO.Info.ExpName;
STs_1unitNstim.Info.Unit=unitINFO.Info.Unit;
STs_1unitNstim.Info.BF_kHz=unitINFO.Info.BF_kHz;
STs_1unitNstim.Info.Threshold_dBSPL=unitINFO.Info.Threshold_dBSPL;
STs_1unitNstim.Info.Q10=unitINFO.Info.Q10;
STs_1unitNstim.Info.SR_sps=unitINFO.Info.SR_sps;
STs_1unitNstim.Info.STMPshift=0;
STs_1unitNstim.STIMtype=STIMtype;

% Make sure there are the proper STIMtype pictures
picList=findPics(STIMtype,[TrackNum,UnitNum]);
if isempty(picList)
    disp(sprintf('*** No %s data for this unit!!!',STIMtype))
    beep
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Get general parameters for this unit, e.g., all BFs, levels, ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load all PICS
PICS=cell(size(picList));
for PICind=1:length(picList)
    PICS{PICind}=loadPic(picList(PICind));
end

if sum(strcmp(STIMtype,{'EHrBFi','EHvNrBFi','TrBFi'}))
    % FIRST, find all octshifts, Features, freqs, polarities,
    % formantsATharmonics
    clear Temp xx
    Temp.octshifts=cell(1,length(picList));
    Temp.FeatureINDs=cell(1,length(picList));
    Temp.freqs_Hz=cell(1,length(picList));
    Temp.FeaturesList=cell(1,length(picList));
    Temp.param_vals=cell(1,length(picList));    % params = levels_dBSPL
    if strcmp(STIMtype,'EHvNrBFi')
        Temp.Vowel_Level_dBSPL=cell(1,length(picList));
    end
    Temp.FirstBADline=Inf+ones(size(picList));
    for PICind=1:length(picList)
        disp(sprintf('      ... Processing picture: %d',picList(PICind)))
        if (PICS{PICind}.General.picture_number~=picList(PICind))|(PICS{PICind}.General.track~=TrackNum)|(PICS{PICind}.General.unit~=UnitNum)
            error(sprintf('data Mismatch for Picture: %d',picList(PICind)));
        end
        
        %% Need to verify no bad-line errors, ow, counting will be off - SO FOR NOW, DISCARD ALL LINES AFTER THIS
        if ~isempty(PICS{PICind}.Stimuli.bad_lines)
            Temp.FirstBADline(PICind)=PICS{PICind}.Stimuli.bad_lines(1);
            disp(sprintf('picture %d: LINE ERRORS, IGNORING ALL LINES STARTING at line #: %d',picList(PICind),Temp.FirstBADline(PICind)));
        end
        
        Temp.octshifts{PICind}=PICS{PICind}.Stimuli.Used.OctShifts_List;
        Temp.FeatureINDs{PICind}=find(strcmp(deblank(PICS{PICind}.Stimuli.Condition.Features{1}),FeaturesText));
        for i=2:length(PICS{PICind}.Stimuli.Condition.Features)
            Temp.FeatureINDs{PICind}=[Temp.FeatureINDs{PICind} find(strcmp(deblank(PICS{PICind}.Stimuli.Condition.Features{i}),FeaturesText))];
        end
        Temp.freqs_Hz{PICind}=PICS{PICind}.Stimuli.Used.FeatureTarget_Hz_List;
        Temp.FeaturesList{PICind}=PICS{PICind}.Stimuli.Used.Features_List;
        %%%% Store params - differs acros stimuli
        if sum(strcmp(STIMtype,{'EHrBFi','TrBFi'}))
            if isfield(PICS{PICind}.Stimuli.Used,'Levels_dBSPL_List')
                Temp.param_vals{PICind}=PICS{PICind}.Stimuli.Used.Levels_dBSPL_List;
            else
                Temp.param_vals{PICind}=PICS{PICind}.Stimuli.Condition.Level_dBSPL;
            end
        elseif strcmp(STIMtype,'EHvNrBFi')
            Temp.param_vals{PICind}=PICS{PICind}.Stimuli.Used.NoiseAttens_dB_List;  % Different for EHvN
            Temp.Vowel_Level_dBSPL{PICind}=PICS{PICind}.Stimuli.Condition.Levels_dBSPL;
        end
        Temp.picNums(PICind)=picList(PICind);
        Temp.HarmIND(PICind)=find(strcmp(deblank(PICS{PICind}.Stimuli.Condition.FormsAtHarmonics),FormsAtHarmonicsText));
        Temp.PolIND(PICind)=find(strcmp(deblank(PICS{PICind}.Stimuli.Condition.InvertPolarity),InvertPolarityText));
    end
    
    TempFeatureINDs=[];
    Tempoctshifts=[];
    Tempparam_vals=[];
    for i=1:length(picList)
        Tempoctshifts=[Tempoctshifts Temp.octshifts{i}];
        TempFeatureINDs=[TempFeatureINDs Temp.FeatureINDs{i}];
        Tempparam_vals=[Tempparam_vals Temp.param_vals{i}];
    end
    if strcmp(STIMtype,'EHvNrBFi')
        TempVowelLevel=[]; % only used for EHvNrBfi
        for i=1:length(picList)
            TempVowelLevel=[TempVowelLevel Temp.Vowel_Level_dBSPL{i}];
        end
        if length(unique(TempVowelLevel))>1
            error('More than 1 vowel level used!!!!')
        end
    end
    SORTocts=unique(Tempoctshifts);
    SORTfeatures=unique(TempFeatureINDs);
    SORTparams=unique(Tempparam_vals);
    SORTpol=unique(Temp.PolIND);
    SORTharms=unique(Temp.HarmIND);
    NumF=length(SORTocts);
    NumP=length(SORTparams);
    NumFEATURES=length(SORTfeatures);
    NumPOL=length(SORTpol);
    NumHARMS=length(SORTharms);
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %% Save parameter values
    STs_1unitNstim.ChannelInfo.channel_parameter='Octave Shift';
    if strcmp(STIMtype,'EHrBFi')
        STs_1unitNstim.ParameterInfo.parameter='Vowel Level';
        STs_1unitNstim.ParameterInfo.parameter_units='dB SPL';
    elseif strcmp(STIMtype,'EHvNrBFi')
        STs_1unitNstim.ParameterInfo.parameter='Noise Attenuation';
        STs_1unitNstim.ParameterInfo.parameter_units='dBatt';
        STs_1unitNstim.ConditionInfo.VowelLevel_dBSPL = TempVowelLevel;
    elseif strcmp(STIMtype,'TrBFi')
        STs_1unitNstim.ParameterInfo.parameter='Tone Level';
        STs_1unitNstim.ParameterInfo.parameter_units='dB SPL';
    end
    STs_1unitNstim.ParameterInfo.param_values=SORTparams;
    STs_1unitNstim.ConditionInfo.features=FeaturesText(SORTfeatures);
    STs_1unitNstim.ConditionInfo.polarities=InvertPolarityText(SORTpol);
    STs_1unitNstim.ConditionInfo.FormsAtHarmonics=FormsAtHarmonicsText(SORTharms);
    STs_1unitNstim.StimInfo.Interleaved='yes';
    for PICind=1:length(picList)
        if PICind==1
            STs_1unitNstim.StimInfo.StimDur_msec=PICS{PICind}.Hardware.Trigger.StmOn;
            STs_1unitNstim.StimInfo.LineDur_msec=PICS{PICind}.Hardware.Trigger.StmOff+PICS{PICind}.Hardware.Trigger.StmOn;
        else
            if STs_1unitNstim.StimInfo.StimDur_msec~=PICS{PICind}.Hardware.Trigger.StmOn;
                error('Mis-matching StimDur');
            end
            if STs_1unitNstim.StimInfo.LineDur_msec~=PICS{PICind}.Hardware.Trigger.StmOn+PICS{PICind}.Hardware.Trigger.StmOff;
                error('Mis-matching LineDur');
            end
        end
    end
    
    STs_1unitNstim.SpikeTrains=cell(NumFEATURES,NumPOL,NumHARMS);
    STs_1unitNstim.ChannelInfo.channel_values=cell(NumFEATURES,NumPOL,NumHARMS);;
    STs_1unitNstim.ChannelInfo.channel_F0_Hz=cell(NumFEATURES,NumPOL,NumHARMS);
    STs_1unitNstim.ChannelInfo.channel_SamplingRate_Hz=cell(NumFEATURES,NumPOL,NumHARMS);
    STs_1unitNstim.ChannelInfo.channel_FeatureTarget_Hz=cell(NumFEATURES,NumPOL,NumHARMS);
    STs_1unitNstim.ChannelInfo.channel_FeatureFreqs_Hz=cell(NumFEATURES,NumPOL,NumHARMS);
    STs_1unitNstim.ChannelInfo.channel_FeatureLevels_dB=cell(NumFEATURES,NumPOL,NumHARMS);
    %%%%%%%%%%%%%%%%%
    %% Find relevant PICnums and excludelines (MAKE SURE TO CHECK unitINFO for
    %% excludelines) and then save SPIKETRAINS
    %%%%%%%%%%%%%%%%%
    
    for FeatIND=SORTfeatures
        disp(sprintf('      ... ... Processing Feature: %s',FeaturesText{FeatIND}))
        for PolIND=1:NumPOL
            for HarmIND=1:NumHARMS
                STs_1unitNstim.ChannelInfo.channel_FeatureFreqs_Hz{find(SORTfeatures==FeatIND),PolIND,HarmIND}=cell(1,NumF);
                STs_1unitNstim.ChannelInfo.channel_FeatureLevels_dB{find(SORTfeatures==FeatIND),PolIND,HarmIND}=cell(1,NumF);
                
                % Find relevant picture indices for each condition (These are Masks [1: yes, 0:no] for each picture in picList)
                clear TempFeatINDs TempHarmINDs TempPolINDs TempINDs
                for i=1:length(picList)
                    TempFeatINDs(i)=sum(ismember(Temp.FeatureINDs{i},FeatIND));
                end
                TempHarmINDs=(Temp.HarmIND==HarmIND);
                TempPolINDs=(Temp.PolIND==PolIND);
                TempINDs=find(TempFeatINDs&TempHarmINDs&TempPolINDs);  % 1 for pictures with current: Feature, Harm, Polarity
                
                if ~isempty(TempINDs)
                    yTEMP.octshifts=SORTocts;
                    yTEMP.freqs_kHz=[];
                    yTEMP.param_vals=SORTparams;
                    
                    %%%%%%%%%%% STORE picNUMS for each Freq and Level
                    STs_1unitNstim.SpikeTrains{find(SORTfeatures==FeatIND),PolIND,HarmIND}=cell(NumF,NumP);
                    for FreqIND=1:NumF
                        disp(sprintf('      ... ... ... Processing Channel: %d',FreqIND))
                        % Store exact freqs for ALL conditions for this Feature (will vary from feature to feature)
                        CONDindFULL=[];
                        for PICind=TempINDs
                            CONDindFULL=[CONDindFULL find((Temp.octshifts{PICind}==SORTocts(FreqIND))&(strcmp(deblank(Temp.FeaturesList{PICind}),FeaturesText{FeatIND})))];
                        end
                        CONDfreq_kHz=unique(Temp.freqs_Hz{PICind}(CONDindFULL)/1000);
                        if length(CONDfreq_kHz)~=1
                            error('Non-unique CONDfreq_kHz, when looking at FreqIND=%d',FreqIND)
                        end
                        yTEMP.freqs_kHz(FreqIND)=CONDfreq_kHz;
                        
                        for ParamIND=1:NumP
                            % Verify unique conditions
                            CONDindFULL=[];
                            for PICind=TempINDs
                                CONDindFULL=[CONDindFULL find((Temp.octshifts{PICind}==SORTocts(FreqIND))&(Temp.param_vals{PICind}==SORTparams(ParamIND))& ...
                                    (strcmp(deblank(Temp.FeaturesList{PICind}),FeaturesText{FeatIND})))];
                            end
                            CONDind=unique(CONDindFULL);
                            if length(CONDind)>1
                                error(sprintf('Non-unique CONDind, when looking at FreqIND=%d and ParamIND=%d',FreqIND,ParamIND))
                            end
                            
                            % Find relevant picture indices for each Freq/Level condition (1: yes, 0:no)
                            clear TempOctINDs TempLevINDs
                            for i=1:length(picList)
                                TempOctINDs(i)=sum(ismember(Temp.octshifts{i},SORTocts(FreqIND)))>0;
                                TempLevINDs(i)=sum(ismember(Temp.param_vals{i},SORTparams(ParamIND)))>0;
                            end
                            TempINDs2=find(TempFeatINDs&TempHarmINDs&TempPolINDs&TempOctINDs&TempLevINDs);  % Mask for Feature,Harm,Polarity,Freq,Level
                            
                            yTEMP.picNums=setdiff(Temp.picNums(TempINDs2),unitINFO.PICS.IgnorePicNums);   % Take out any IgnorePICS from unitINFO
                            % Store excludeLines for each picture, based on CONDind and TOTALconds
                            yTEMP.excludeLines=cell(size(TempINDs2));
                            ii=0;
                            for i=TempINDs2
                                ii=ii+1;
                                TOTALconds=length(PICS{i}.Stimuli.Used.FeatureTarget_Hz_List);
                                GOODlines=CONDind:TOTALconds:PICS{i}.Stimuli.fully_presented_lines;
                                %% Take out any lines beyond badlines
                                GOODlines=GOODlines(find(GOODlines<Temp.FirstBADline(i)));
                                IGNORElines=setdiff(1:PICS{i}.Stimuli.fully_presented_lines,GOODlines);
                                
                                % Check for BADlines from unitINFO
                                currPIC = Temp.picNums(i);
                                IGNORElines = union(IGNORElines,unitINFO.PICS.picBADlines{find(unitINFO.PICS.PIClist==currPIC)});
                                
                                yTEMP.excludeLines{ii}=IGNORElines;
                            end
                            
                            % Save spikeTrain for this Freq/Level
                            PIC=concatPICS_STMP(yTEMP.picNums,yTEMP.excludeLines);
                            STs_1unitNstim.SpikeTrains{find(SORTfeatures==FeatIND),PolIND,HarmIND}{FreqIND,ParamIND}=PIC.x.spikes{1};
                            
                            % Store Fundamental Frequency and SamplingRate for each channel
                            if ParamIND==1
                                STs_1unitNstim.ChannelInfo.channel_values{find(SORTfeatures==FeatIND),PolIND,HarmIND}=SORTocts;
                                STs_1unitNstim.ChannelInfo.channel_F0_Hz{find(SORTfeatures==FeatIND),PolIND,HarmIND}(FreqIND)=PIC.FundamentalFreq_Hz;
                                STs_1unitNstim.ChannelInfo.channel_SamplingRate_Hz{find(SORTfeatures==FeatIND),PolIND,HarmIND}(FreqIND)=PIC.x.Stimuli.Used.UpdateRate_Hz_List(CONDind);
                                STs_1unitNstim.ChannelInfo.channel_FeatureTarget_Hz{find(SORTfeatures==FeatIND),PolIND,HarmIND}(FreqIND)=PIC.x.Stimuli.Used.FeatureTarget_Hz_List(CONDind);
                                STs_1unitNstim.ChannelInfo.channel_FeatureFreqs_Hz{find(SORTfeatures==FeatIND),PolIND,HarmIND}{FreqIND}=PIC.FeatureFreqs_Hz;
                                if strcmp(STIMtype,'TrBFi')
                                    STs_1unitNstim.ChannelInfo.channel_FeatureLevels_dB{find(SORTfeatures==FeatIND),PolIND,HarmIND}{FreqIND}=0;
                                else
                                    STs_1unitNstim.ChannelInfo.channel_FeatureLevels_dB{find(SORTfeatures==FeatIND),PolIND,HarmIND}{FreqIND}=PIC.FeatureLevels_dB;
                                end
                            else
                                if STs_1unitNstim.ChannelInfo.channel_F0_Hz{find(SORTfeatures==FeatIND),PolIND,HarmIND}(FreqIND)~=PIC.FundamentalFreq_Hz
                                    error('Fundamental Freqs dont match')
                                end
                                if STs_1unitNstim.ChannelInfo.channel_SamplingRate_Hz{find(SORTfeatures==FeatIND),PolIND,HarmIND}(FreqIND)~=PIC.x.Stimuli.Used.UpdateRate_Hz_List(CONDind)
                                    error('Sampling Rates dont match')
                                end
                                if sum(STs_1unitNstim.ChannelInfo.channel_FeatureFreqs_Hz{find(SORTfeatures==FeatIND),PolIND,HarmIND}{FreqIND}~=PIC.FeatureFreqs_Hz)
                                    error('Feature Freqs dont match')
                                end
                                if ~strcmp(STIMtype,'TrBFi')
                                    if sum(STs_1unitNstim.ChannelInfo.channel_FeatureLevels_dB{find(SORTfeatures==FeatIND),PolIND,HarmIND}{FreqIND}~=PIC.FeatureLevels_dB)
                                        error('Feature Freqs dont match')
                                    end
                                end
                                
                            end
                            
                        end  % end Param
                    end  % end Freq
                end
            end % End Invert Polarity
        end % End FormsatHarms
    end % End: FeatIND
end % End: STIMtype


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADDED June 5, 2009 - Boley
if sum(strcmp(STIMtype,{'BBN_SX','BBNrBFi','WAVreBFi'}))
    BASEsamprate_Hz=33003.300330;
    % FIRST, find all stim, Fs, levels, and polarities USED
    clear Temp xx
    Temp.StimNames=cell(1,length(picList));
    Temp.SampRates=cell(1,length(picList));
    Temp.Attens=cell(1,length(picList));
    Temp.FirstBADline=Inf+ones(size(picList));
    for PICind=1:length(picList)
        disp(sprintf('      ... Processing picture: %d',picList(PICind)))
        if (PICS{PICind}.General.picture_number~=picList(PICind))|(PICS{PICind}.General.track~=TrackNum)|(PICS{PICind}.General.unit~=UnitNum)
            error(sprintf('data Mismatch for Picture: %d',picList(PICind)));
        end
        
        %% Need to verify no bad-line errors, ow, counting will be off - SO FOR NOW, DISCARD ALL LINES AFTER THIS
        if ~isempty(PICS{PICind}.Stimuli.bad_lines)
            Temp.FirstBADline(PICind)=PICS{PICind}.Stimuli.bad_lines(1);
            disp(sprintf('picture %d: LINE ERRORS, IGNORING ALL LINES STARTING at line #: %d',picList(PICind),Temp.FirstBADline(PICind)));
        end
        
        Temp.StimNames{PICind}=PICS{PICind}.Stimuli.list;
        for j=1:length(Temp.StimNames{PICind})
            [pathstr,Temp.StimNames{PICind}{j},ext,versn]=fileparts(Temp.StimNames{PICind}{j});
            %% Set based on specific file names used for BBN_SX = BBN_A,
            %% BBN_AN, BBN_B, BBN_BN - REMOVE negative polarity versions
            if strcmp(Temp.StimNames{PICind}{j}(end),'N')
                Temp.StimNames{PICind}{j}=Temp.StimNames{PICind}{j}(1:end-1);
            end
        end
        Temp.SampRates{PICind}=PICS{PICind}.Stimuli.updateRate_Hz;
        Temp.Attens{PICind}=PICS{PICind}.Stimuli.attens;
        Temp.picNums(PICind)=picList(PICind);
    end
    
    TempStimNames={};
    TempSampRates=[];
    TempAttens=[];
    for i=1:length(picList)
        for j=1:length(Temp.StimNames{i})
            TempStimNames{end+1}=Temp.StimNames{i}{j};
        end
        TempSampRates=[TempSampRates Temp.SampRates{i}];
        TempAttens=[TempAttens Temp.Attens{i}];
    end
    SORTstim=unique(TempStimNames);
    SORTsamprates=unique(TempSampRates);
    SORTocts=log2(SORTsamprates./BASEsamprate_Hz);
    SORTattens=unique(TempAttens);
    NumSTIM=length(SORTstim);
    NumSAMPRATE=length(SORTsamprates);
    NumATTEN=length(SORTattens);
    SORTpols=InvertPolarityText;
    NumPOL=2;
end

%%%%%%%%%%%%%%%%%%%%%%%%
%% Save parameter values
if sum(strcmp(STIMtype,{'BBN_SX','BBNrBFi','WAVreBFi'}))
    
    %% %%%%%%%%%%%%%%%%%%%
    % Updated 1/14/09: M. Heinz for WAVreBFi
    % This book-keeping assumes:
    %   - can be different WAV files (all have +&- versions, with the
    %     neg-polarity filename with an "...N" at the end
    %   - ALL wav files use the same Oct-shifts
    %   - ***Each wav file may have a different set of Attenuations, thus some
    %     attenuations in list may be empty for a given STIM
    %     *** This is the big change from BBNrBFi, but now works for both
    
    STs_1unitNstim.ChannelInfo.channel_parameter='Octave Shift';
    STs_1unitNstim.ParameterInfo.parameter='Attenuation';
    STs_1unitNstim.ParameterInfo.parameter_units='dBatt';
    STs_1unitNstim.ParameterInfo.param_values=SORTattens;
    STs_1unitNstim.ConditionInfo.features=SORTstim;  % Here, features are the different stimuli
    STs_1unitNstim.ConditionInfo.polarities=InvertPolarityText;
    STs_1unitNstim.StimInfo.Interleaved='yes';
    for PICind=1:length(picList)
        if PICind==1
            STs_1unitNstim.StimInfo.StimDur_msec=PICS{PICind}.Hardware.Trigger.StmOn;
            STs_1unitNstim.StimInfo.LineDur_msec=PICS{PICind}.Hardware.Trigger.StmOff+PICS{PICind}.Hardware.Trigger.StmOn;
        else
            if STs_1unitNstim.StimInfo.StimDur_msec~=PICS{PICind}.Hardware.Trigger.StmOn;
                error('Mis-matching StimDur');
            end
            if STs_1unitNstim.StimInfo.LineDur_msec~=PICS{PICind}.Hardware.Trigger.StmOn+PICS{PICind}.Hardware.Trigger.StmOff;
                error('Mis-matching LineDur');
            end
        end
    end
    STs_1unitNstim.SpikeTrains=cell(NumSTIM,NumPOL);
    STs_1unitNstim.ChannelInfo.channel_values=cell(NumSTIM,NumPOL);
    STs_1unitNstim.ChannelInfo.channel_SamplingRate_Hz=cell(NumSTIM,NumPOL);
    %%%%%%%%%%%%%%%%%
    %% Find relevant PICnums and excludelines (MAKE SURE TO CHECK unitINFO for
    %% excludelines) and then save SPIKETRAINS
    %%%%%%%%%%%%%%%%%
    
    for StimIND=1:NumSTIM
        for PolIND=1:NumPOL
            %% Need to adjust condition FILENAME back to include 'N' at end it
            %% needed - for BBN_SX
            condFILENAME=SORTstim{StimIND};
            if strcmp(InvertPolarityText{PolIND},'yes')
                condFILENAME(end+1)='N';
            end
            disp(sprintf('      ... ... Processing Stimuli: %s, PolarityInversion: %s',condFILENAME,upper(InvertPolarityText{PolIND})))
            
            % Find relevant picture indices for each condition (These are Masks [1: yes, 0:no] for each picture in picList)
            clear TempStimINDs TempHarmINDs TempPolINDs TempINDs
            for i=1:length(picList)
                TempStimINDs(i)=sum(strcmp(Temp.StimNames{i},SORTstim{StimIND}))>0;
            end
            TempPolINDs=ones(size(TempStimINDs)); % by definition, all BBN_SX PICs have +&-
            TempINDs=find(TempStimINDs&TempPolINDs);  % 1 for pictures with current: Feature, Polarity
            
            if ~isempty(TempINDs)
                %%%%%%%%%%% STORE picNUMS for each Freq and Level
                STs_1unitNstim.SpikeTrains{StimIND,PolIND}=cell(NumSAMPRATE,NumATTEN);
                for srIND=1:NumSAMPRATE
                    clear TempOctINDs
                    disp(sprintf('      ... ... ... Processing Sampling Rate: %.2f Hz [Octave Shift = %.2f]',SORTsamprates(srIND),SORTocts(srIND)))
                    % Find relevant picture indices for each SampRate condition (1: yes, 0:no)
                    for i=1:length(picList)
                        TempOctINDs(i)=sum(ismember(Temp.SampRates{i},SORTsamprates(srIND)))>0;
                    end
                    
                    for ParamIND=1:NumATTEN
                        disp(sprintf('      ... ... ... ... Processing Attenuation : %.2f dB',SORTattens(ParamIND)))
                        % Find relevant picture indices for each Level condition (1: yes, 0:no)
                        clear TempLevINDs
                        for i=1:length(picList)
                            % 							TempLevINDs(i)=sum(ismember(Temp.Attens{i},SORTattens(ParamIND)))>0;
                            % Need Atten to exist for this stim specifically!  %
                            % 1/14/09
                            TempLevINDs(i)=sum((ismember(Temp.Attens{i},SORTattens(ParamIND)))&(strcmp(Temp.StimNames{i},SORTstim{StimIND})));  %%
                            
                            %% TO DO LATER - CLEANUP - why does this have to be so
                            %% complicated, one param at a time?? Why can't we just
                            %% inside the final for loop do one TempINDs for all 4
                            %% params???
                        end
                        TempINDs2=find(TempStimINDs&TempPolINDs&TempOctINDs&TempLevINDs);  % Mask for Feature,Polarity,Freq,Level
                        
                        picNums=setdiff(Temp.picNums(TempINDs2),unitINFO.PICS.IgnorePicNums);   % Take out any IgnorePICS from unitINFO
                        TempINDs2=find(ismember(Temp.picNums,picNums));	%re-adjust TempINDs2
                        
                        if ~isempty(TempINDs2)
                            % Store excludeLines for each picture, based on CONDind and TOTALconds
                            excludeLines=cell(size(TempINDs2));
                            %% Find relevant lines numbers for each PICTURE for this
                            %% condition, and exclude the rest
                            ii=0;
                            for i=TempINDs2
                                IGNORElines=[];
                                ii=ii+1;
                                
                                
                                %%%% REWRITE as CONDind=N, then find reps of N -
                                %%%% see STMPsaveSTs file
                                
                                % Find current filename
                                for j=1:length(PICS{i}.Line.file)
                                    [pathstr,tempNAME,ext,versn]=fileparts(PICS{i}.Line.file{j});
                                    if ~strcmp(tempNAME,condFILENAME)
                                        IGNORElines=[IGNORElines j];
                                    elseif (j>=Temp.FirstBADline(i))|(j>PICS{i}.Stimuli.fully_presented_lines)   %% Take out any lines beyond badlines or fully presented lines
                                        IGNORElines=[IGNORElines j];
                                    end
                                end
                                % Find current samp rate
                                SRlistTEMP=repmat(PICS{i}.Stimuli.Used.UpdateRate_Hz_List,1,PICS{i}.Stimuli.repetitions);
                                SRlistTEMP=SRlistTEMP(1:PICS{i}.spikes{1}(end,1));
                                IGNORElines2=find(SRlistTEMP~=SORTsamprates(srIND));
                                IGNORElines=union(IGNORElines,IGNORElines2);
                                
                                % Check for BADlines from unitINFO
                                currPIC = Temp.picNums(i);
                                IGNORElines = union(IGNORElines,unitINFO.PICS.picBADlines{find(unitINFO.PICS.PIClist==currPIC)});
                                excludeLines{ii}=IGNORElines;
                                disp(sprintf('      ... ... ... ... ... Taking spikes from PIC: %d [TOTAL REPS included = %d; 1st included line = %d]', ...
                                    currPIC,length(PICS{i}.Line.file)-length(excludeLines{ii}),min(setdiff(1:length(PICS{i}.Line.file),excludeLines{ii}))))
                            end
                            
                            % Save spikeTrain for this Freq/Level
                            PIC=concatPICS_BBNSAC(picNums,excludeLines);
                            STs_1unitNstim.SpikeTrains{StimIND,PolIND}{srIND,ParamIND}=PIC.x.spikes{1};
                        else
                            disp('      ... ... ... ... ... ***** NO SPIKES for this condition *****')
                            STs_1unitNstim.SpikeTrains{StimIND,PolIND}{srIND,ParamIND}=[];
                        end
                        
                        % Store Fundamental Frequency and SamplingRate for each channel
                        if ParamIND==1
                            STs_1unitNstim.ChannelInfo.channel_values{StimIND,PolIND}=SORTocts;
                            STs_1unitNstim.ChannelInfo.channel_SamplingRate_Hz{StimIND,PolIND}=SORTsamprates;
                        else
                            if ~isempty(TempINDs2)
                                if STs_1unitNstim.ChannelInfo.channel_SamplingRate_Hz{StimIND,PolIND}(srIND)~=PIC.x.Stimuli.updateRate_Hz
                                    error('Sampling Rates dont match')
                                end
                            end
                        end
                        
                    end  % end Param
                end  % end SampRate
            end % end isempty
        end % End Invert Polarity
    end % End: StimIND
end % END:STIMtype
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SAVE SpikeTrains
disp(['   *** Saving new file: "' STs_1unitNstim_filename '" ...'])
eval(['save ''' fullfile(STMPanal_dir,STs_1unitNstim_filename) ''' STs_1unitNstim'])

return;
