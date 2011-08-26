function [SACSCCfunctions,SACSCCmetrics,paramsOUT] = CCCanal(SpikeTrains,paramsIN,PLOT_15panel)
% File: CCCanal.m
% M. Heinz/ J. Swaminathan
% May 20, 2008
%
% General function for neural cross-correlation analyses/metrics.
%
% this program computes the envelope and fine structure coding and
% correlations from SAC/SCC/XpAC/XpCC analyses of neurophysiological spike
% train data.
%
% SpikeTrains: cell array with 4 spiketrains {condition (1,2), polarity (plus,minus)}
% paramsIN(OUT): defaults used if not included
%         .SACbinwidth_msec
%         .ignoreONSET_msec
%         .dur1_msec
%         .dur2_msec
% PLOT_15panel: 1=yes, plot fill 15-panel figure; 0=don't plot
%
% SACSCCfunctions: structure with all functions returned
%                .delays_usec
%                .SAC_A_plus   : SAC for condition A+
%                .SAC_A_minus  : SAC for condition A-
%                .SAC_A_avg    : SAC for AVG(A+,A-)
%                .SAC_B_plus   : SAC for condition B+
%                .SAC_B_minus  : SAC for condition B-
%                .SAC_B_avg    : SAC for AVG(B+,B-)
%                .SCC_AB_plus  : SCC for condition A+/B+
%                .SCC_AB_minus : SCC for condition A-/B-
%                .SCC_AB_avg   : SCC for AVG(A+/B+,A-/B-)
%                .XpAC_A_plus  : XpAC for condition A+/A-
%                .XpAC_A_minus : XpAC for condition A-/A+
%                .XpAC_A_avg   : XpAC for AVG(A+/A-,A-/A+)
%                .XpAC_B_plus  : XpAC for condition B+/B-
%                .XpAC_B_minus : XpAC for condition B-/B+
%                .XpAC_B_avg   : XpAC for AVG(B+/B-,B-/B+)
%                .XpCC_AB_plus : XpCC for condition A+/B-
%                .XpCC_AB_minus: XpCC for condition A-/B+
%                .XpCC_AB_avg  : XpCC for AVG(A+/B-,A-/B+)
%                .SUMCOR_A     : SUMCOR for condition A = AVG(SAC_A_avg,XpAC_A_avg)
%                .SUMCOR_B     : SUMCOR for condition B = AVG(SAC_B_avg,XpAC_B_avg)
%                .SUMCOR_AB    : SUMCOR for condition 3 = AVG(SCC_A_avg,XpCC_A_avg)
%                .DIFCOR_A     : DIFCOR for condition A = SAC_A_avg-XpAC_A_avg
%                .DIFCOR_B     : DIFCOR for condition B = SAC_B_avg-XpAC_B_avg
%                .DIFCOR_AB    : DIFCOR for condition AB = SCC_AB_avg-XpCC_AB_avg
%
% SACSCCmetrics: structure with all summary metrics returned
%                .CCCenv1
%                .CCCtfs
%                .CDenv_usec
%                .CDtfs_usec
%                .CDscc_usec
%                .DCpeak
%                .SCpeak_1
%                .NumDrivenSpikes(2,2)
%                .AvgRate_sps(2,2)

LOADexistingDATA=0;

if ~LOADexistingDATA

	if ~exist('PLOT_15panel','var'),    PLOT_15panel=1;   end

	%% set the parameters for SAC and SCC %%%%%%%
	paramsOUT=paramsIN;
	% ignore ONSET - set default if not specified
	if ~isfield(paramsOUT,'SCC_onsetIGNORE_sec')
		paramsOUT.SCC_onsetIGNORE_sec=0.05;
	end
	% smoothing for functions
	TriFiltWidthSAC=5; TriFiltWidthDC=5; TriFiltWidthSC=1;
	paramsOUT.TriFiltWidthSAC=TriFiltWidthSAC;
	paramsOUT.TriFiltWidthDC=TriFiltWidthDC;
	paramsOUT.TriFiltWidthSC=TriFiltWidthSC;
	% SCC binwidth (Joris, 2003)
	paramsOUT.DELAYbinwidth_sec=50e-6;
	% FFT window length - to get 1-Hz sampling
	Nfft_psd=round(1/paramsOUT.DELAYbinwidth_sec);  % 20000
	paramsOUT.Nfft_psd=Nfft_psd;
	% SAC/SCC window length
	MAXdelay_sec=0.0125;              paramsOUT.MAXdelay_sec=MAXdelay_sec;
	% stim duration for SAC analysis = min of both conditions
	stimdur_sec=min(paramsIN.durA_msec,paramsIN.durB_msec)/1000;
	paramsOUT.stimdur_sec=stimdur_sec;
	% Limit the number of spikes in the SAC/SCC analysis
	paramsOUT.MAXspikes=5000;  % If used, windowSTs will cut out extra REPs beyond MAXspikes - used to avoid memory limits
	% minimum DCpeak for A and B to compute CCCtfs
	paramsOUT.minDCpeak_CCCtfs=0.1;
	
	%% Window spike times (ignore 1st 50 ms, and cut to shorter of 2 stim)
	[ST_A_plus,NumDrivenSpikes(1,1)]=windowSTs(SpikeTrains{1,1},paramsOUT.SCC_onsetIGNORE_sec,stimdur_sec,paramsOUT.MAXspikes);
	[ST_A_minus,NumDrivenSpikes(1,2)]=windowSTs(SpikeTrains{1,2},paramsOUT.SCC_onsetIGNORE_sec,stimdur_sec,paramsOUT.MAXspikes);
	[ST_B_plus,NumDrivenSpikes(2,1)]=windowSTs(SpikeTrains{2,1},paramsOUT.SCC_onsetIGNORE_sec,stimdur_sec,paramsOUT.MAXspikes);
	[ST_B_minus,NumDrivenSpikes(2,2)]=windowSTs(SpikeTrains{2,2},paramsOUT.SCC_onsetIGNORE_sec,stimdur_sec,paramsOUT.MAXspikes);
	paramsOUT.SCCdur_sec=stimdur_sec-paramsOUT.SCC_onsetIGNORE_sec;
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp(sprintf('NUMBER OF SPIKES (A_plus, A_minus, B_plus, B_minus) = %d, %d, %d, %d',NumDrivenSpikes(1,1),NumDrivenSpikes(1,2), ...
		NumDrivenSpikes(2,1),NumDrivenSpikes(2,2)))
	disp(sprintf('NUMBER OF REPS (A_plus, A_minus, B_plus, B_minus) = %d, %d, %d, %d',length(ST_A_plus),length(ST_A_minus), ...
		length(ST_B_plus),length(ST_B_minus)))

	disp('   ... Computing SACs/SCCs for AN spikes ...')
	[SACSCCs,AvgRate_sps]=SACSCCanal(ST_A_plus,ST_A_minus,ST_B_plus,ST_B_minus,paramsOUT);

	if ~isnan(NumDrivenSpikes(2,1))
		
			%% RANDOMIZE spike times within same window used above (ignore 1st 50 ms, and cut to shorter of 2 stim)
			STrand_A_plus=randomizeSTs(ST_A_plus,paramsOUT.SCC_onsetIGNORE_sec,stimdur_sec);
			STrand_A_minus=randomizeSTs(ST_A_minus,paramsOUT.SCC_onsetIGNORE_sec,stimdur_sec);
			STrand_B_plus=randomizeSTs(ST_B_plus,paramsOUT.SCC_onsetIGNORE_sec,stimdur_sec);
			STrand_B_minus=randomizeSTs(ST_B_minus,paramsOUT.SCC_onsetIGNORE_sec,stimdur_sec);

			disp(sprintf('   ... Computing SACs/SCCs for RANDOM spikes (for bias removal) ...'))
			[SACSCCs_rand,AvgRate_rand_sps]=SACSCCanal(STrand_A_plus,STrand_A_minus,STrand_B_plus,STrand_B_minus,paramsOUT);
			
	else
		SACSCCs_rand=[];
	end

% 	save testSACSCC.mat
else
	load testSACSCC.mat
	paramsOUT.CF_Hz=1.0*1000;
	
end

%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TODO 5/20/08
% *0.5) SETUP BBN_A, BBN_B vs DURATION
% *1) clean up function, verify returning all - fix returning bookkeeping
% *2) ADD WINDOW correcton in function to all ROW1 functions
% *3) setup random spikes, call same function, discard most, except PSD,CSD
% for RS
%
%
%% 5/22/08 - 11pm
% %1) clean up plotting, book-keeping
% *- add all metrics to structure, to allow outside plotting function for
% *15-panel (SACSCCfunctions, SACSCCmetrics)
%
% *4) finish CD metrics, add to plot
% *verify with some examples!!! (picking of peaks)
%
% 4.75) verify FIG 2 in JON paper
%  * - go over with Ganesh - FRIDAY
%
%% 6) MOVE ON to DELTA CF - VERIFY STMP FOR EXP!!!!!!!
% 6.5) set up to run only COL 1 for [CFi, CFi]
% 7) clean up code for LAB use






%% PSD/CSD summations - Compute summed energy in ENVELOPE spectral
%% densities over various frequency ranges
% find INDs in freqVEC_Hz for relevant cutoffs
paramsOUT.CCCenv_LOWmod_Hz=10;
paramsOUT.CCCenv_HIGHmod_Hz=300;  %?? min(CF,300)???  TRY SMALLER WINDOW 100 Hz resolut?
[y,CCCenv_ZEROmod_ind]=min(abs(SACSCCs.freqVEC-0));
[y,CCCenv_LOWmod_ind]=min(abs(SACSCCs.freqVEC-paramsOUT.CCCenv_LOWmod_Hz));
[y,CCCenv_HIGHmod_ind]=min(abs(SACSCCs.freqVEC-paramsOUT.CCCenv_HIGHmod_Hz));
[y,CCCenv_CFmod_ind]=min(abs(SACSCCs.freqVEC-paramsOUT.CF_Hz));
[y,CCCenv_NYQmod_ind]=min(abs(SACSCCs.freqVEC-0.5*(1/paramsOUT.DELAYbinwidth_sec)));

% CCCenv based frequency range (10-300 Hz)
SACSCCmetrics.sums.sumPSDsc_A_CCCenv = sum(SACSCCs.PSDsc_A(CCCenv_LOWmod_ind:CCCenv_HIGHmod_ind));
if ~isnan(NumDrivenSpikes(2,1))
	SACSCCmetrics.sums.sumPSDsc_B_CCCenv = sum(SACSCCs.PSDsc_B(CCCenv_LOWmod_ind:CCCenv_HIGHmod_ind));
	SACSCCmetrics.sums.sumCSDsc_AB_CCCenv = sum(SACSCCs.CSDsc_AB(CCCenv_LOWmod_ind:CCCenv_HIGHmod_ind));
	SACSCCmetrics.sums.sumPSDsc_Arand_CCCenv = sum(SACSCCs_rand.PSDsc_A(CCCenv_LOWmod_ind:CCCenv_HIGHmod_ind));
	SACSCCmetrics.sums.sumPSDsc_Brand_CCCenv = sum(SACSCCs_rand.PSDsc_B(CCCenv_LOWmod_ind:CCCenv_HIGHmod_ind));
	SACSCCmetrics.sums.sumCSDsc_ABrand_CCCenv = sum(SACSCCs_rand.CSDsc_AB(CCCenv_LOWmod_ind:CCCenv_HIGHmod_ind));
end

% 0-CCCenv_HIGH based frequency range (0-300 Hz: for SCpeak3)
SACSCCmetrics.sums.sumPSDsc_A_0_CCCenv = sum(SACSCCs.PSDsc_A(CCCenv_ZEROmod_ind:CCCenv_HIGHmod_ind));
if ~isnan(NumDrivenSpikes(2,1))
	SACSCCmetrics.sums.sumPSDsc_B_0_CCCenv = sum(SACSCCs.PSDsc_B(CCCenv_ZEROmod_ind:CCCenv_HIGHmod_ind));
	SACSCCmetrics.sums.sumCSDsc_AB_0_CCCenv = sum(SACSCCs.CSDsc_AB(CCCenv_ZEROmod_ind:CCCenv_HIGHmod_ind));
	SACSCCmetrics.sums.sumPSDsc_Arand_0_CCCenv = sum(SACSCCs_rand.PSDsc_A(CCCenv_ZEROmod_ind:CCCenv_HIGHmod_ind));
	SACSCCmetrics.sums.sumPSDsc_Brand_0_CCCenv = sum(SACSCCs_rand.PSDsc_B(CCCenv_ZEROmod_ind:CCCenv_HIGHmod_ind));
	SACSCCmetrics.sums.sumCSDsc_ABrand_0_CCCenv = sum(SACSCCs_rand.CSDsc_AB(CCCenv_ZEROmod_ind:CCCenv_HIGHmod_ind));
end

% 0-CF frequency range
SACSCCmetrics.sums.sumPSDsc_A_0_CF = sum(SACSCCs.PSDsc_A(CCCenv_ZEROmod_ind:CCCenv_CFmod_ind));
if ~isnan(NumDrivenSpikes(2,1))
	SACSCCmetrics.sums.sumPSDsc_B_0_CF = sum(SACSCCs.PSDsc_B(CCCenv_ZEROmod_ind:CCCenv_CFmod_ind));
	SACSCCmetrics.sums.sumCSDsc_AB_0_CF = sum(SACSCCs.CSDsc_AB(CCCenv_ZEROmod_ind:CCCenv_CFmod_ind));
	SACSCCmetrics.sums.sumPSDsc_Arand_0_CF = sum(SACSCCs_rand.PSDsc_A(CCCenv_ZEROmod_ind:CCCenv_CFmod_ind));
	SACSCCmetrics.sums.sumPSDsc_Brand_0_CF = sum(SACSCCs_rand.PSDsc_B(CCCenv_ZEROmod_ind:CCCenv_CFmod_ind));
	SACSCCmetrics.sums.sumCSDsc_ABrand_0_CF = sum(SACSCCs_rand.CSDsc_AB(CCCenv_ZEROmod_ind:CCCenv_CFmod_ind));
end

% CCCenvLOW-CF frequency range (10-CF: for CCCenv4)
SACSCCmetrics.sums.sumPSDsc_A_CCCenv_CF = sum(SACSCCs.PSDsc_A(CCCenv_LOWmod_ind:CCCenv_CFmod_ind));
if ~isnan(NumDrivenSpikes(2,1))
	SACSCCmetrics.sums.sumPSDsc_B_CCCenv_CF = sum(SACSCCs.PSDsc_B(CCCenv_LOWmod_ind:CCCenv_CFmod_ind));
	SACSCCmetrics.sums.sumCSDsc_AB_CCCenv_CF = sum(SACSCCs.CSDsc_AB(CCCenv_LOWmod_ind:CCCenv_CFmod_ind));
	SACSCCmetrics.sums.sumPSDsc_Arand_CCCenv_CF = sum(SACSCCs_rand.PSDsc_A(CCCenv_LOWmod_ind:CCCenv_CFmod_ind));
	SACSCCmetrics.sums.sumPSDsc_Brand_CCCenv_CF = sum(SACSCCs_rand.PSDsc_B(CCCenv_LOWmod_ind:CCCenv_CFmod_ind));
	SACSCCmetrics.sums.sumCSDsc_ABrand_CCCenv_CF = sum(SACSCCs_rand.CSDsc_AB(CCCenv_LOWmod_ind:CCCenv_CFmod_ind));
end

% 0-Fs/2 frequency range
SACSCCmetrics.sums.sumPSDsc_A_0_NYQ = sum(SACSCCs.PSDsc_A(CCCenv_ZEROmod_ind:CCCenv_NYQmod_ind));
if ~isnan(NumDrivenSpikes(2,1))
	SACSCCmetrics.sums.sumPSDsc_B_0_NYQ = sum(SACSCCs.PSDsc_B(CCCenv_ZEROmod_ind:CCCenv_NYQmod_ind));
	SACSCCmetrics.sums.sumCSDsc_AB_0_NYQ = sum(SACSCCs.CSDsc_AB(CCCenv_ZEROmod_ind:CCCenv_NYQmod_ind));
	SACSCCmetrics.sums.sumPSDsc_Arand_0_NYQ = sum(SACSCCs_rand.PSDsc_A(CCCenv_ZEROmod_ind:CCCenv_NYQmod_ind));
	SACSCCmetrics.sums.sumPSDsc_Brand_0_NYQ = sum(SACSCCs_rand.PSDsc_B(CCCenv_ZEROmod_ind:CCCenv_NYQmod_ind));
	SACSCCmetrics.sums.sumCSDsc_ABrand_0_NYQ = sum(SACSCCs_rand.CSDsc_AB(CCCenv_ZEROmod_ind:CCCenv_NYQmod_ind));
end

%% Compute metrics
%% Peak Heights
% SAC peaks - this is CI from Joris et al 2006 (HR) [don't subtract 1]
SACSCCmetrics.SACpeak_A=max(SACSCCs.SAC_A_avg);
if ~isnan(NumDrivenSpikes(2,1))
	SACSCCmetrics.SACpeak_B=max(SACSCCs.SAC_B_avg);
	SACSCCmetrics.SCCpeak_AB=max(SACSCCs.SCC_AB_avg);
end
% DIFCOR peaks
SACSCCmetrics.DCpeak_A=max(SACSCCs.DIFCOR_A);
if ~isnan(NumDrivenSpikes(2,1))
	SACSCCmetrics.DCpeak_B=max(SACSCCs.DIFCOR_B);
	SACSCCmetrics.DCpeak_AB=max(SACSCCs.DIFCOR_AB);
end
% SUMCOR peaks - don't subtract 1 (Louage et al 2004)
% SCpeak2: raw peaks
SACSCCmetrics.SCpeak2_A=max(SACSCCs.SUMCOR_A);
if ~isnan(NumDrivenSpikes(2,1))
	SACSCCmetrics.SCpeak2_B=max(SACSCCs.SUMCOR_B);
	SACSCCmetrics.SCpeak2_AB=max(SACSCCs.SUMCOR_AB);
end
%% USE THIS ONE: SCpeak: Eq. 2 = adjusted SC peak height [0,CF]/[0,Fs/2]
SACSCCmetrics.SCpeak_A=1+(SACSCCmetrics.SCpeak2_A-1)*SACSCCmetrics.sums.sumPSDsc_A_0_CF/SACSCCmetrics.sums.sumPSDsc_A_0_NYQ;
if ~isnan(NumDrivenSpikes(2,1))
	SACSCCmetrics.SCpeak_B=1+(SACSCCmetrics.SCpeak2_B-1)*SACSCCmetrics.sums.sumPSDsc_B_0_CF/SACSCCmetrics.sums.sumPSDsc_B_0_NYQ;
	SACSCCmetrics.SCpeak_AB=1+(SACSCCmetrics.SCpeak2_AB-1)*SACSCCmetrics.sums.sumCSDsc_AB_0_CF/SACSCCmetrics.sums.sumCSDsc_AB_0_NYQ;
end
% SCpeak_3: alternative adjusted SC peak height [0,300]/[0,Fs/2]
SACSCCmetrics.SCpeak3_A=1+(SACSCCmetrics.SCpeak2_A-1)*SACSCCmetrics.sums.sumPSDsc_A_0_CCCenv/SACSCCmetrics.sums.sumPSDsc_A_0_NYQ;
if ~isnan(NumDrivenSpikes(2,1))
	SACSCCmetrics.SCpeak3_B=1+(SACSCCmetrics.SCpeak2_B-1)*SACSCCmetrics.sums.sumPSDsc_B_0_CCCenv/SACSCCmetrics.sums.sumPSDsc_B_0_NYQ;
	SACSCCmetrics.SCpeak3_AB=1+(SACSCCmetrics.SCpeak2_AB-1)*SACSCCmetrics.sums.sumCSDsc_AB_0_CCCenv/SACSCCmetrics.sums.sumCSDsc_AB_0_NYQ;
end

if ~isnan(NumDrivenSpikes(2,1))
	
	%% Neural Cross Correlation Coefficients
	if (SACSCCmetrics.DCpeak_A>=paramsOUT.minDCpeak_CCCtfs)&(SACSCCmetrics.DCpeak_B>=paramsOUT.minDCpeak_CCCtfs)
		SACSCCmetrics.CCCtfs = SACSCCmetrics.DCpeak_AB/(sqrt(SACSCCmetrics.DCpeak_A*SACSCCmetrics.DCpeak_B));
	else
		SACSCCmetrics.CCCtfs = NaN;
	end
	%% USE THIS ONE: Eq. 3 - 10-300 PSD, with noise bias removed
	SACSCCmetrics.CCCenv=max([0 (SACSCCmetrics.sums.sumCSDsc_AB_CCCenv-SACSCCmetrics.sums.sumCSDsc_ABrand_CCCenv)])/ ...
		sqrt((SACSCCmetrics.sums.sumPSDsc_A_CCCenv-SACSCCmetrics.sums.sumPSDsc_Arand_CCCenv)* ...
		(SACSCCmetrics.sums.sumPSDsc_B_CCCenv-SACSCCmetrics.sums.sumPSDsc_Brand_CCCenv));
	% raw SC peaks
	SACSCCmetrics.CCCenv2 = (SACSCCmetrics.SCpeak2_AB-1)/(sqrt((SACSCCmetrics.SCpeak2_A-1)*(SACSCCmetrics.SCpeak2_B-1)));
	% 10-300 Hz, bias NOT removed
	SACSCCmetrics.CCCenv3=SACSCCmetrics.sums.sumCSDsc_AB_CCCenv/sqrt(SACSCCmetrics.sums.sumPSDsc_A_CCCenv*SACSCCmetrics.sums.sumPSDsc_B_CCCenv);
	% 10-CF Hz, noise bias removed
	SACSCCmetrics.CCCenv4=max([0 (SACSCCmetrics.sums.sumCSDsc_AB_CCCenv_CF-SACSCCmetrics.sums.sumCSDsc_ABrand_CCCenv_CF)])/ ...
		sqrt((SACSCCmetrics.sums.sumPSDsc_A_CCCenv_CF-SACSCCmetrics.sums.sumPSDsc_Arand_CCCenv_CF)* ...
		(SACSCCmetrics.sums.sumPSDsc_B_CCCenv_CF-SACSCCmetrics.sums.sumPSDsc_Brand_CCCenv_CF));
	% adjusted SC peaks
	SACSCCmetrics.CCCenv5 = (SACSCCmetrics.SCpeak_AB-1)/(sqrt((SACSCCmetrics.SCpeak_A-1)*(SACSCCmetrics.SCpeak_B-1)));

	%% Characteristic Delays
	%% NEED TO IMPROVE - find largest, unless there is a "close-to-largest" that
	% is closer to zero!
	% probably use some criterion for 2nd largest peak (if with 5%)???
	SACSCCmetrics.CDscc_usec=findCD_SCC(SACSCCs.SCC_AB_avg,SACSCCs.delays_usec);
	SACSCCmetrics.CDenv_usec=findCD_SCC(SACSCCs.SUMCOR_AB,SACSCCs.delays_usec);
	SACSCCmetrics.CDtfs_usec=findCD_SCC(SACSCCs.DIFCOR_AB,SACSCCs.delays_usec);
end

%% Book-keeping
SACSCCfunctions=SACSCCs;
if ~isnan(NumDrivenSpikes(2,1))
	SACSCCfunctions.rand.PSDsc_A=SACSCCs_rand.PSDsc_A;
	SACSCCfunctions.rand.PSDsc_B=SACSCCs_rand.PSDsc_B;
	SACSCCfunctions.rand.CSDsc_AB=SACSCCs_rand.CSDsc_AB;
end

SACSCCmetrics.NumDrivenSpikes=NumDrivenSpikes;
SACSCCmetrics.AvgRate_sps=AvgRate_sps;

if PLOT_15panel
% 	paramsOUT.XLIMIT_delay=5;
% 	paramsOUT.XLIMIT_PSDhigh=700;
% 	paramsOUT.YLIMIT_SClow=.5;
	
	plot_CCCanal(SACSCCfunctions,SACSCCmetrics,paramsOUT)
end



return;

function [SACSCCs,AvgRate_sps]=SACSCCanal(ST_A_plus,ST_A_minus,ST_B_plus,ST_B_minus,paramsOUT);

MAXdelay_ind=round(paramsOUT.MAXdelay_sec/paramsOUT.DELAYbinwidth_sec);  % old XLIM=250

%% COLUMN 1: CONDITION=A
%% Compute SAC (A+ and A-) %%%%%%%%%%%%%%%%%%%%
[SAC_A_plus,delays_usec,AvgRate_sps(1,1),xxNspikes] = ShufAutoCorr(ST_A_plus,paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
[SAC_A_minus,delays_usec,AvgRate_sps(1,2),xxNspikes] = ShufAutoCorr(ST_A_minus,paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
% SAME for all conditions
ZEROind=find(delays_usec==0); 
SACinds=(ZEROind-MAXdelay_ind:ZEROind+MAXdelay_ind);
SACdelays_usec=delays_usec(SACinds);
% SAC limited-data Window CORRECTION 
TEMP=linspace(1,0,ZEROind);
WindowCORRECTION=[TEMP fliplr(TEMP(1:end-1))];
% % USED TO TAKE OUT WINDOW CORRECTION to test
% beep
% disp('NO WINDOW CORRECTION APPLIED')
% WindowCORRECTION=zeros(size(WindowCORRECTION));

SAC_A_plus=SAC_A_plus+WindowCORRECTION;
SAC_A_minus=SAC_A_plus+WindowCORRECTION;
% smooth and AVG correlograms
SAC_A_plus=trifilt(SAC_A_plus(SACinds),paramsOUT.TriFiltWidthSAC);
SAC_A_minus=trifilt(SAC_A_minus(SACinds),paramsOUT.TriFiltWidthSAC);
SAC_A_avg=(SAC_A_plus+SAC_A_minus)/2;

%% Compute XpAC (A+/A- and A-/A+) %%%%%%%%%%%%%%%%%%%%
[XpAC_A_plus,delays_usec,xxAVGrate,xxNspikes] = ShufCrossCorr({ST_A_plus ST_A_minus},paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
% SAC limited-data Window CORRECTION 
XpAC_A_plus=XpAC_A_plus+WindowCORRECTION;
% smooth and AVG correlograms
XpAC_A_plus=trifilt(XpAC_A_plus(SACinds),paramsOUT.TriFiltWidthSAC);
XpAC_A_minus=fliplr(XpAC_A_plus);
XpAC_A_avg=(XpAC_A_plus+XpAC_A_minus)/2;

%% Compute DIFCOR and SUMCOR for Condition A 
DIFCOR_A=trifilt(SAC_A_avg-XpAC_A_avg,paramsOUT.TriFiltWidthDC);  % Diffcorr =  Avg(SAC) - Avg(XpAC)
SUMCOR_A=trifilt((SAC_A_avg+XpAC_A_avg)/2,paramsOUT.TriFiltWidthSC); %Sumcorr = 1/2(Avg(SAC) + Avg(XpAC))

%% Compute ENVELOPE POWER SPECTRAL DENSITY for Condition A
PSDsc_A=abs(fft((SUMCOR_A-1),paramsOUT.Nfft_psd));
freqVEC=(0:length(PSDsc_A)-1)/length(PSDsc_A)*paramsOUT.Nfft_psd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(ST_B_plus)
	%% COLUMN 2: CONDITION=B
	%% Compute SAC (B+ and B-) %%%%%%%%%%%%%%%%%%%%
	[SAC_B_plus,delays_usec,AvgRate_sps(2,1),xxNspikes] = ShufAutoCorr(ST_B_plus,paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
	[SAC_B_minus,delays_usec,AvgRate_sps(2,2),xxNspikes] = ShufAutoCorr(ST_B_minus,paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
	% SAC limited-data Window CORRECTION
	SAC_B_plus=SAC_B_plus+WindowCORRECTION;
	SAC_B_minus=SAC_B_minus+WindowCORRECTION;
	% smooth and AVG correlograms
	SAC_B_plus=trifilt(SAC_B_plus(SACinds),paramsOUT.TriFiltWidthSAC);
	SAC_B_minus=trifilt(SAC_B_minus(SACinds),paramsOUT.TriFiltWidthSAC);
	SAC_B_avg=(SAC_B_plus+SAC_B_minus)/2;

	%% Compute XpAC (B+/B- and B-/B+) %%%%%%%%%%%%%%%%%%%%
	[XpAC_B_plus,delays_usec,xxAVGrate,xxNspikes] = ShufCrossCorr({ST_B_plus ST_B_minus},paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
	% SAC limited-data Window CORRECTION
	XpAC_B_plus=XpAC_B_plus+WindowCORRECTION;
	% smooth and AVG correlograms
	XpAC_B_plus=trifilt(XpAC_B_plus(SACinds),paramsOUT.TriFiltWidthSAC);
	XpAC_B_minus=fliplr(XpAC_B_plus);
	XpAC_B_avg=(XpAC_B_plus+XpAC_B_minus)/2;

	%% Compute DIFCOR and SUMCOR for Condition B
	DIFCOR_B=trifilt(SAC_B_avg-XpAC_B_avg,paramsOUT.TriFiltWidthDC);  % Diffcorr =  Avg(SAC) - Avg(XpAC)
	SUMCOR_B=trifilt((SAC_B_avg+XpAC_B_avg)/2,paramsOUT.TriFiltWidthSC); %Sumcorr = 1/2(Avg(SAC) + Avg(XpAC))

	%% Compute ENVELOPE POWER SPECTRAL DENSITY for Condition B
	PSDsc_B=abs(fft((SUMCOR_B-1),paramsOUT.Nfft_psd));
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	%% COLUMN 3: CONDITION=AB
	%% Compute SCC (A+/B+ and A-/B-) %%%%%%%%%%%%%%%%%%%%
	[SCC_AB_plus,delays_usec,xxAVGrate,xxNspikes] = ShufCrossCorr({ST_A_plus ST_B_plus},paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
	[SCC_AB_minus,delays_usec,xxAVGrate,xxNspikes] = ShufCrossCorr({ST_A_minus ST_B_minus},paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
	% SAC limited-data Window CORRECTION
	SCC_AB_plus=SCC_AB_plus+WindowCORRECTION;
	SCC_AB_minus=SCC_AB_minus+WindowCORRECTION;
	% smooth and AVG correlograms
	SCC_AB_plus=trifilt(SCC_AB_plus(SACinds),paramsOUT.TriFiltWidthSAC);
	SCC_AB_minus=trifilt(SCC_AB_minus(SACinds),paramsOUT.TriFiltWidthSAC);
	SCC_AB_avg=(SCC_AB_plus+SCC_AB_minus)/2;

	%% Compute XpCC (A+/B- and A-/B+) %%%%%%%%%%%%%%%%%%%%
	[XpCC_AB_plus,delays_usec,xxAVGrate,xxNspikes] = ShufCrossCorr({ST_A_plus ST_B_minus},paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
	[XpCC_AB_minus,delays_usec,xxAVGrate,xxNspikes] = ShufCrossCorr({ST_A_minus ST_B_plus},paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
	% SAC limited-data Window CORRECTION
	XpCC_AB_plus=XpCC_AB_plus+WindowCORRECTION;
	XpCC_AB_minus=XpCC_AB_minus+WindowCORRECTION;
	% smooth and AVG correlograms
	XpCC_AB_plus=trifilt(XpCC_AB_plus(SACinds),paramsOUT.TriFiltWidthSAC);
	XpCC_AB_minus=trifilt(XpCC_AB_minus(SACinds),paramsOUT.TriFiltWidthSAC);
	XpCC_AB_avg=(XpCC_AB_plus+XpCC_AB_minus)/2;

	%% Compute DIFCOR and SUMCOR for Condition AB
	DIFCOR_AB=trifilt(SCC_AB_avg-XpCC_AB_avg,paramsOUT.TriFiltWidthDC);  % Diffcorr =  Avg(SCC) - Avg(XpCC)
	SUMCOR_AB=trifilt((SCC_AB_avg+XpCC_AB_avg)/2,paramsOUT.TriFiltWidthSC); %Sumcorr = 1/2(Avg(SCC) + Avg(XpCC))

	%% Compute ENVELOPE CROSS SPECTRAL DENSITY for Condition AB
	CSDsc_AB=abs(fft((SUMCOR_AB-1),paramsOUT.Nfft_psd));
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
	SAC_B_plus=[];
	AvgRate_sps(2,1)=NaN;
	SAC_B_minus=[];
	AvgRate_sps(2,2)=NaN;
	SAC_B_avg=[];
	XpAC_B_plus=[];
	XpAC_B_minus=[];
	XpAC_B_avg=[];
	DIFCOR_B=[];
	SUMCOR_B=[];
	PSDsc_B=[];
	SCC_AB_plus=[];
	SCC_AB_minus=[];
	SCC_AB_avg=[];
	XpCC_AB_plus=[];
	XpCC_AB_minus=[];
	XpCC_AB_avg=[];
	DIFCOR_AB=[];
	SUMCOR_AB=[];
	CSDsc_AB=[];
end

%% Store all params, output for return;
% % % SACSCCfunctions: structure with ALL (including plus and minus) functions returned
% % SACSCCs=struct('delays_usec',SACdelays_usec,'SAC_A_plus',SAC_A_plus,'SAC_A_minus',SAC_A_minus,'SAC_A_avg',SAC_A_avg,'SAC_B_plus',SAC_B_plus, ...
% % 	'SAC_B_minus',SAC_B_minus,'SAC_B_avg',SAC_B_avg,'SCC_AB_plus',SCC_AB_plus,'SCC_AB_minus',SCC_AB_minus,'SCC_AB_avg',SCC_AB_avg, ...
% % 	'XpAC_A_plus',XpAC_A_plus,'XpAC_A_minus',XpAC_A_minus,'XpAC_A_avg',XpAC_A_avg,'XpAC_B_plus',XpAC_B_plus,'XpAC_B_minus',XpAC_B_minus, ...
% % 	'XpAC_B_avg',XpAC_B_avg,'XpCC_AB_plus',XpCC_AB_plus,'XpCC_AB_minus',XpCC_AB_minus,'XpCC_AB_avg',XpCC_AB_avg,'SUMCOR_A',SUMCOR_A, ...
% % 	'SUMCOR_B',SUMCOR_B,'SUMCOR_AB',SUMCOR_AB,'DIFCOR_A',DIFCOR_A,'DIFCOR_B',DIFCOR_B,'DIFCOR_AB',DIFCOR_AB,'PSDsc_A',PSDsc_A,'PSDsc_B',PSDsc_B, ...
% % 	'CSDsc_AB',CSDsc_AB,'freqVEC',freqVEC);
% SACSCCfunctions: structure with ALL useful functions returned
SACSCCs=struct('delays_usec',SACdelays_usec,'SAC_A_avg',SAC_A_avg,'SAC_B_avg',SAC_B_avg,'SCC_AB_avg',SCC_AB_avg,'XpAC_A_avg',XpAC_A_avg, ...
	'XpAC_B_avg',XpAC_B_avg,'XpCC_AB_avg',XpCC_AB_avg,'SUMCOR_A',SUMCOR_A, ...
	'SUMCOR_B',SUMCOR_B,'SUMCOR_AB',SUMCOR_AB,'DIFCOR_A',DIFCOR_A,'DIFCOR_B',DIFCOR_B,'DIFCOR_AB',DIFCOR_AB,'PSDsc_A',PSDsc_A,'PSDsc_B',PSDsc_B, ...
	'CSDsc_AB',CSDsc_AB,'freqVEC',freqVEC);

return;

