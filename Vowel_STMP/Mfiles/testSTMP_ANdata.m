% File: testSTMP_ANdata.m
% M. Heinz Jun 08, 2007
% Analyzes one unit all the way through the STMP process.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global STMP_dir STMP_ExpList
% global FeaturesText FormsAtHarmonicsText InvertPolarityText


% close previous RLF plots
OPENfigs=get(0,'Children');
close(intersect([14 1 2 3 4 100 101],OPENfigs))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Verify parameters and experiment, unit are valid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear ExpDate
TESTunitNUM=7;  % 1: ARO (111804, 1.28), 2: ISH (041805, 2.04), 3: Purdue1 (111606, 2.09); 4: Purdue2 (041307, 2.03), 
% 99: ISH 2009 (122208, 6.17)
%%%% Specify ExpDate if not provided
%%%
%% CHECKED: 20Jun2007  (UpTO RPS)
% TESTunitNUM: *1, *2, *3, *4, *5
if ~exist('ExpDate','var')
   if TESTunitNUM==1
      ExpDate='111804'; UnitName='1.28'; STIMtype='EHrBFi'; DATAtype='Nchans1stim';  % ARO2005: EHrBFI, with EHrlv
   elseif TESTunitNUM==2
      ExpDate='041805'; UnitName='2.04'; STIMtype='EHvNrBFi'; DATAtype='Nchans1stim'; % ISH2006: EHvnrBFI, with EHrlv and EHINrlv
	elseif TESTunitNUM==3
		ExpDate='111606'; UnitName='2.09'; STIMtype='EHvNrBFi'; DATAtype='Nchans1stim';  % Purdue #1
	elseif TESTunitNUM==4
		ExpDate='041307'; UnitName='2.03'; STIMtype='EHvNrBFi'; DATAtype='Nchans1stim';  % Purdue #2

	elseif TESTunitNUM==5
      ExpDate='111804'; UnitName='1.29'; STIMtype='EHrBFi'; DATAtype='Nchans1stim';  % ARO2005: EHrBFI, with EHrlv
	elseif TESTunitNUM==6
      ExpDate='121404'; UnitName='1.03'; STIMtype='EHrBFi'; DATAtype='Nchans1stim';  % ARO2005: EHrBFI, with EHrlv
	elseif TESTunitNUM==7
      ExpDate='121404'; UnitName='1.07'; STIMtype='EHrBFi'; DATAtype='Nchans1stim';  % ARO2005: EHrBFI, with EHrlv
	elseif TESTunitNUM==8
      ExpDate='121404'; UnitName='2.04'; STIMtype='EHrBFi'; DATAtype='Nchans1stim';  % ARO2005: EHrBFI, with EHrlv
	elseif TESTunitNUM==9
      ExpDate='121404'; UnitName='2.09'; STIMtype='EHrBFi'; DATAtype='Nchans1stim';  % ARO2005: EHrBFI, with EHrlv
   elseif TESTunitNUM==10
      ExpDate='041805'; UnitName='2.04'; STIMtype='TrBFi'; DATAtype='Nchans1stim';  % TONE re BFI, with Tone RLV
      
   elseif TESTunitNUM==99
      ExpDate='122208'; UnitName='6.17'; STIMtype='WAVreBFi'; DATAtype='Nchans1stim';  % Broadband Noise
   end            
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
disp(sprintf('\n*********** STMP: Processing STIMtype: %s for:  Experiment: ''%s''; Unit: %d.%02d  **********',STIMtype,ExpName,TrackNum,UnitNum))

% STMPsave_UnitINFO(TrackNum,UnitNum,1)
% Make sure GEN works here too, try to make this a GENERAL function!
% 6/17/08 - CCC analyis for ASA
GENsave_UnitINFO(TrackNum,UnitNum,1)

STMPsave_STs_1unitNstim_interleaved(ExpDate,UnitName,STIMtype)

STMPconvert_STs_Nchans1stim(ExpDate,UnitName,STIMtype)

STMPplot_RLFs(ExpDate,UnitName)

% STMPcompute_PERhists(ExpDate,UnitName,STIMtype,DATAtype)
% 
% STMPplot_PERhists(ExpDate,UnitName,STIMtype,DATAtype)
% 
% STMPcompute_DFTs_RSPALSR(ExpDate,UnitName,STIMtype,DATAtype)
% STMPplot_DFTs(ExpDate,UnitName,STIMtype,DATAtype)
% STMPplot_RPSs(ExpDate,UnitName,STIMtype,DATAtype)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CHECKED: 20Jun2007  (UpTO RPS)
% TESTunitNUM: *1, *2, *3, *4, *5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% NEEDS to be verified (compare w/ ISH 2007):
STMPplot_SMPs(ExpDate,UnitName,STIMtype,DATAtype)

%% NOT fully tested yet:
% STMPcompute_SCCs(ExpDate,UnitName,STIMtype,DATAtype);
try STMPBBNcompute_SCCs(ExpDate,UnitName,STIMtype,DATAtype,'BBN_A_',20), end
pack

%% Needs to be written (see ISH 2009)
% STMPplot_SCCs(ExpDate,UnitName,STIMtype,DATAtype)





%%%%%%%%% 
% *1) compute and save PerHists (separate file)
% *2) plot PerHists (as in ARO - with general function)
% *3) plot rate, synch, phase
% 4) compute and save PSTs (separate file)
% *5) compute and save DFTs (separate file)
% 6) SETUP REST 
%    - ISH
%    - Purdue
%    - Tone reBF


%%%%%%%%%%%%%%%
%%% 6/8/07 TODO
% *1) Setup unitINFO with JUST basic info to save
%       - SR, Q10, BF, PICs, BADlines (run this FOR EACH UNIT at set up - DONT SKIP!)
%       see NOTES, see makeNOHRdatalist - setup to just runn for each unit
%    GET RID OF WHOLE EXP DataList - too confusing - just do 1 unit at a
%    time!
% 2) Setup STsaving (ST: spike train)
%   *a) orig
%   *b) shifted
% *3) Setup  - PerHist computing, plotting 
%%     General plot for NEUROGRAMS
% *4) Setup  - DFT computing, plotting (& R,S,P,ALSR)
% 3) Setup  - PSTs (for SACs) computing, plotting 
% 3) Setup  - SCC,SAC computing, plotting 



