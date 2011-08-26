function ReviewUnitTriggering(TrackNum,UnitNum,checkRASTER)
% function ReviewUnitTriggering(Track,Unit)
% Created: M. Heinz 30Dec2005
% For: R03
%
% Goes through all pictures for a unit and reviews: Triggering, Errors, and
% comments

if ~exist('checkRASTER')
	checkRASTER=1;
end
	
picList=findPics('*',[TrackNum,UnitNum]);

for PICnum=picList
	x=loadPic(PICnum);
	disp(sprintf('Picture #: %d, filename: %s',PICnum,getfileName(PICnum)))

	if isfield(x.General,'trigger')
		disp(sprintf('   Trigger: %s',upper(x.General.trigger)))
		if sum(strcmp(deblank(x.General.trigger),{'Poor','Fair'}))
			beep
		end
	end
	
	if isfield(x.General,'comment')
		disp(sprintf('   Comment: %s',upper(x.General.comment)))
	end
	
	if isfield(x.General,'run_errors')
		for i=1:length(x.General.run_errors)
			if ~sum(strcmp(x.General.run_errors{i}, ...
					{'In function ''DALinloop_NI_wavfiles'': Input waveform ', ...
					'has been truncated to fit requested duration. ', ...
					'has been repeated to fill requested duration. '}))
			disp(sprintf('   Run_errors: %s',x.General.run_errors{i}))
			end 
		end
	end

	if checkRASTER
		if ~strcmp('tc',getTAG(getfileName(PICnum)))
			PICview(PICnum,'')
			[SR_sps,SRpics,SRests_sps,SRdata] = PIC_calcSR(PICnum);
			disp(sprintf('MEAN SPONT RATE = %.1f sp/sec',SR_sps))
			beep
			input('Press Enter to move to next PICTURE');

		
		end
	end	
	
end

	

return;



