function PIC=concatPICS_BBNSAC(picNums,excludeLines)
% File: PIC=concatPICS_BBNSAC(picNums,excludeLines)
% Jun 17, 2008
% From: PIC=concatPICS_STMP(picNums,excludeLines)
% M.Heinz 10Sep2004.  Taken from PSTview.m
% Concatenates multiple pst-style pics; allows specification of lines to be excluded for each picture.
% Assumes current directory is the data directory
% Error checking is based on NOHR stimulus conditions
% Also, stores Feature Freqs, BF, and Fixed Frequency for NOHR stimulus conditions
%
% Updated 07Jan2005: to allow for global storage of PICS
%
% Usage: PIC=concatPICS_NOHR(picNums,excludeLines)
% picNums: vector of picture numbers
% excludeLines: cell array with vectors for each picture with any lines to be excluded
%
%% TO DO - clean up all commented lines taken out 6/18/08 from STMP 
%

global SavedPICS SavedPICnums SavedPICSuse

NUMpics=length(picNums);
if ~exist('excludeLines','var')
   excludeLines=cell(1,NUMpics);
end

% If empty, it means we're not using this feature
if isempty(SavedPICSuse)
   SavedPICSuse=0;
end
% If we're not using this feature, clear any left over data
if ~SavedPICSuse
   SavedPICS=[];
   SavedPICnums=[];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load all pictures into PICS
PICS=cell(1,NUMpics);
for PICind=1:NUMpics
   PICS{PICind}.num = picNums(PICind);
   PICS{PICind}.excludeLines = excludeLines{PICind};
   
   % If we're not saving PICS, just load PIC as usual
   if ~SavedPICSuse
      [PICS{PICind}.x,errorMSG]=loadPic(PICS{PICind}.num);
      if ~isempty(errorMSG)
         error(errorMSG);
      end
   else % We're saving PICS
      if ~sum(ismember(SavedPICnums,PICS{PICind}.num))
         [PICS{PICind}.x,errorMSG]=loadPic(PICS{PICind}.num);
         if ~isempty(errorMSG)
            error(errorMSG);
         end
         % SAVE PIC in global memory if SavedPICS is being used         
         SavedPICS{length(SavedPICS)+1}=PICS{PICind}.x;
         SavedPICnums(length(SavedPICnums)+1)=PICS{PICind}.num;
      else % we have this PIC stored already
         PICS{PICind}.x=SavedPICS{find(ismember(SavedPICnums,PICS{PICind}.num))};
      end      
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup PIC with all data, and individual run_errors, triggers, comments
PIC=PICS{1};  % Take all data from PICS{1}, and check rest against this
PIC.TEMPLATE=getTAG(getFileName(PIC.num));
if strcmp(PIC.TEMPLATE,'SX')
	PIC.TEMPLATE='BBN_SX';  % to fix the mistake of using underscotre in TEMPLATE name: BBN_SX!!
end
PIC=rmfield(PIC,'excludeLines');
PIC.nums=zeros(1,length(PICS));
PIC.excludeLines=cell(1,length(PICS));
PIC.triggers=cell(1,length(PICS));
PIC.comments=cell(1,length(PICS));
PIC.runerrors=cell(1,length(PICS));
for PICind=1:length(PICS)
   PIC.nums(PICind)=PICS{PICind}.num;
   PIC.excludeLines{PICind}=PICS{PICind}.excludeLines;
   PIC.triggers{PICind}=PICS{PICind}.x.General.trigger;
   PIC.comments{PICind}=PICS{PICind}.x.General.comment;
   PIC.run_errors{PICind}=PICS{PICind}.x.General.run_errors;
end
PIC=rmfield(PIC,'num');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Verify UNIT AND conditions are the same for all PICS 
for PICind=2:length(PICS)
	% Verify TEMPLATE is the same
	tempTEMPLATE=getTAG(getFileName(PIC.nums(PICind)));
	if strcmp(tempTEMPLATE,'SX')
		tempTEMPLATE='BBN_SX';  % to fix the mistake of using underscotre in TEMPLATE name: BBN_SX!!
	end
	if ~strcmp(PIC.TEMPLATE,tempTEMPLATE)
		error(sprintf('Template mismatch between pictures %d (%s) and %d (%s)',PIC.nums(1),PIC.TEMPLATE, ...
			PIC.nums(PICind),tempTEMPLATE))
	end
   
   % Verify Unit # is the same
   if ((PIC.x.General.track~=PICS{PICind}.x.General.track) | (PIC.x.General.unit~=PICS{PICind}.x.General.unit))
      error(sprintf('Unit # mismatch between pictures %d and %d',PIC.nums(1),PIC.nums(PICind)))      
	end

	% Verify Conditions are the same
%    if strcmp(PIC.TEMPLATE,'TrBF')
%       if ((PIC.x.Stimuli.Condition.BaseFrequency_kHz~=PICS{PICind}.x.Stimuli.Condition.BaseFrequency_kHz) | ...
%             (PIC.x.Stimuli.Condition.FreqOffset_octs~=PICS{PICind}.x.Stimuli.Condition.FreqOffset_octs) | ...
%             ~strcmp(PIC.x.Stimuli.Condition.Offset_Direction,PICS{PICind}.x.Stimuli.Condition.Offset_Direction) | ...
%             (PIC.x.Stimuli.Condition.Level_dBSPL~=PICS{PICind}.x.Stimuli.Condition.Level_dBSPL))
%          error(sprintf('TrBF Condition mismatch between pictures %d and %d',PIC.nums(1),PIC.nums(PICind)))
%       end
%    elseif strcmp(PIC.TEMPLATE,'TrFF')
%       if ((PIC.x.Stimuli.Condition.BaseFrequency_kHz~=PICS{PICind}.x.Stimuli.Condition.BaseFrequency_kHz) | ...
%             (PIC.x.Stimuli.Condition.Level_dBSPL~=PICS{PICind}.x.Stimuli.Condition.Level_dBSPL))
%          error(sprintf('TrFF Condition mismatch between pictures %d and %d',PIC.nums(1),PIC.nums(PICind)))
%       end
%    elseif strcmp(PIC.TEMPLATE,'EHrBF')
%       if ((PIC.x.Stimuli.Condition.BaseFrequency_kHz~=PICS{PICind}.x.Stimuli.Condition.BaseFrequency_kHz) | ...
%             (PIC.x.Stimuli.Condition.FreqOffset_octs~=PICS{PICind}.x.Stimuli.Condition.FreqOffset_octs) | ...
%             ~strcmp(PIC.x.Stimuli.Condition.Offset_Direction,PICS{PICind}.x.Stimuli.Condition.Offset_Direction) | ...
%             ~strcmp(PIC.x.Stimuli.Condition.Feature,PICS{PICind}.x.Stimuli.Condition.Feature) | ...
%             ~strcmp(PIC.x.Stimuli.Condition.FormsAtHarmonics,PICS{PICind}.x.Stimuli.Condition.FormsAtHarmonics) | ...
%             ~strcmp(PIC.x.Stimuli.Condition.InvertPolarity,PICS{PICind}.x.Stimuli.Condition.InvertPolarity) | ...
%             (PIC.x.Stimuli.Condition.Level_dBSPL~=PICS{PICind}.x.Stimuli.Condition.Level_dBSPL))
%          error(sprintf('EHrBF Condition mismatch between pictures %d and %d',PIC.nums(1),PIC.nums(PICind)))
%       end
%    elseif strcmp(PIC.TEMPLATE,'EHrFF')
%       if ((PIC.x.Stimuli.Condition.BaseFrequency_kHz~=PICS{PICind}.x.Stimuli.Condition.BaseFrequency_kHz) | ...
%             ~strcmp(PIC.x.Stimuli.Condition.Feature,PICS{PICind}.x.Stimuli.Condition.Feature) | ...
%             ~strcmp(PIC.x.Stimuli.Condition.FormsAtHarmonics,PICS{PICind}.x.Stimuli.Condition.FormsAtHarmonics) | ...
%             ~strcmp(PIC.x.Stimuli.Condition.InvertPolarity,PICS{PICind}.x.Stimuli.Condition.InvertPolarity) | ...
%             (PIC.x.Stimuli.Condition.Level_dBSPL~=PICS{PICind}.x.Stimuli.Condition.Level_dBSPL))
%          error(sprintf('EHrFF Condition mismatch between pictures %d and %d',PIC.nums(1),PIC.nums(PICind)))
%       end
%    elseif sum(strcmp(PIC.TEMPLATE,{'EHrBFi','EHvNrBFi','TrBFi','SACrlv'}))
%       disp(sprintf('NEED TO UPDATE concatPICS_STMP - check conditions are the same for TEMPLATE: %s',PIC.TEMPLATE))
% 	end
	disp(sprintf('NEED TO UPDATE concatPICS_BBNSAC - check conditions are the same for TEMPLATE: %s',PIC.TEMPLATE))
end

%%%%%%%%%%%%%%%%%%%%%%%%
% Concatenate all spikes
NumLines=0;
PIC.x.spikes{1}=[];
for PICind=1:length(PICS)
   x=PICS{PICind}.x;
   % verify all spikes in PICS{1} are from fully presented lines
   LASTgoodSpike=max(find(x.spikes{1}(:,1)==x.Stimuli.fully_presented_lines));
   if isempty(LASTgoodSpike)  % if last spike is an earlier line than all presented
      LASTgoodSpike=length(x.spikes{1}(:,1));
   end
   SpikeINDs=find(~ismember(x.spikes{1}(1:LASTgoodSpike,1)',PICS{PICind}.excludeLines));
   % Adjust Line Numbers, if needed
   x.spikes{1}(1:end,1)=x.spikes{1}(1:end,1)+NumLines;
   PIC.x.spikes{1}=[PIC.x.spikes{1};x.spikes{1}(SpikeINDs,:)];
   NumLines=PIC.x.spikes{1}(end,1);
end
%%%%%%% Take out gaps in line numbers
PIC.x.spikes{1}(:,1)=cumsum([1 sign(diff(PIC.x.spikes{1}(:,1)))'])';
NumLines=PIC.x.spikes{1}(end,1);
PIC.x.Stimuli.fully_presented_lines=NumLines;  % update for concatenated PIC

return;
