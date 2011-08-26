function [SR_sps,SRpics,SRests_sps,SRdata] = PIC_calcSR(PICnum)
% File: PIC_calcSR.m

% TrackNum=str2num(unitINFO.Info.Unit(1));
% UnitNum=str2num(unitINFO.Info.Unit(3:4));

SRpics=[];
SRests_sps=[];
SRdata={};

% UnitPicList=findPics('*',[TrackNum,UnitNum]);
% TCpicList=findPics('tc',[TrackNum,UnitNum]);


[SR_sps,lineSRs_sps]=calcSRpic(PICnum);  % No excludelines for Now

%% LATER: setup multiple PIC computations

if ~isnan(SR_sps)
   SRpics=[SRpics PICnum];
   SRests_sps=[SRests_sps SR_sps];
   SRdata{length(SRdata)+1}=lineSRs_sps;
end

