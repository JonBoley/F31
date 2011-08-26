function [SACSCCfunctions,SACSCCmetrics,paramsOUT] = CCCanal_3(SpikeTrains,paramsIN,PLOT_15panel)
% File: CCCanal.m
% M. Heinz/ J. Swaminathan
% May 20, 2008
%
% Modified: July 29, 2008 to generalize SCadjpeak and CCCenv comps based on
% different frequency regions.
% *** MODIFIED CCCanal_3 ONLY, 0, 1, 2 have not been updated with:
%           1) max([0 NFremoval]) to imaginary numbers from sqrt(<0)
%           2) generalization to compute several frequency ranges
%
% CCCanal_0: 1 basic RUN using all spike REPs
% CCCanal_1: Use all spike REPs, run 5 RS reps, avg PSDs, compute 1 CCC 
% CCCanal_2: Use all spike REPs, run 5 indep RUNS, avg 5 CCCs 
% CCCanal_3: BOOTSTRAP 5 RUNS using 80% reps each, run RS rep for each,
%        3A: compute individual CCCs for each RUN, then AVG to one CCC
%        3B: AVG 5 PSDs (AN and RS), the compute 1 CCC (store as RUN: 6)
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

if ~exist('PLOT_15panel','var'),    PLOT_15panel=1;   end

%% set the parameters for SAC and SCC %%%%%%%
paramsOUT=paramsIN;
% ignore ONSET - set default if not specified
if ~isfield(paramsOUT,'SCC_onsetIGNORE_sec')
	paramsOUT.SCC_onsetIGNORE_sec=0.05;
end
% Default SCpeak and CCCenvs to use
if ~isfield(paramsOUT,'SCpeak_TOUSE')
	paramsOUT.SCpeak_TOUSE='IFFTraw';
end
if ~isfield(paramsOUT,'CCCenv_TOUSE')
	paramsOUT.CCCenv_TOUSE='0-300, subBIAS';
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
if ~isfield(paramsOUT,'MAXspikes')
	paramsOUT.MAXspikes=3700;  % If used, windowSTs will cut out extra REPs beyond MAXspikes - used to avoid memory limits
end
% minimum DCpeak for A and B to compute CCCtfs
paramsOUT.minDCpeak_CCCtfs=0.1;
% CF: use minimum of 2 CFs to be most conservative for adjSCpeak
CF_Hz=min([paramsOUT.CF_A_Hz paramsOUT.CF_B_Hz]);  
paramsOUT.SACSCC_CF_Hz=CF_Hz;
%% BOOTSTRAPPING params
paramsOUT.BOOTSTRAP_percentdata=0.80;  % Percent of AN reps to include in each RUN

%% BOOTSTRAPPING - SETUP spike REPS to use for each RUN
for i=1:2
	for j=1:2
		Nreps(i,j)=length(SpikeTrains{i,j});
		Ninclude(i,j)=ceil(paramsOUT.BOOTSTRAP_percentdata*Nreps(i,j)); % some rounding issues with ceil/floor - have to do it this way
		Nexclude(i,j)=max([1 Nreps(i,j)-Ninclude(i,j)]);  % at least 1 rep excluded to make sure runs are different
		Navgs(i,j)=floor(Nreps(i,j)/Nexclude(i,j));  % Number of RUNS for Averaging CCCs, or PSDs
	end
end
paramsOUT.BOOTSTRAP_Navgs=min(min(Navgs));  % Number of RUNS for Averaging CCCs, or PSDs
for i=1:2
	for j=1:2
		BOOTinds{i,j}=cell(1,paramsOUT.BOOTSTRAP_Navgs);
	end
end

for BOOTrep1=1:paramsOUT.BOOTSTRAP_Navgs
	for i=1:2
		for j=1:2
			TEMPinds=setdiff((1:Nreps(i,j)),(Nreps(i,j)-Nexclude(i,j)+1:Nreps(i,j))-(BOOTrep1-1)*Nexclude(i,j));
			%% NEED TO RANDOMIZE SPIKES to avoid repllication if > 5000 spikes
			%% in 50
			%% reps (e.g., if only 20 reps are needed to get 5000 spikes, then
			%% 1st 3 bootstrap samples will all use 1:20
			randINDs=randperm(length(TEMPinds));
			BOOTinds{i,j}{BOOTrep1}=TEMPinds(randINDs);
		end
	end
end


%% First paramsOUT.BOOTSTRAP_Navgs reps are boot strap reps using 80% of
%% reps, last BOOTrep uses the AVG of all reps for all computations
for BOOTrep=1:paramsOUT.BOOTSTRAP_Navgs+1
	% Run BOOTSTRAP_Navgs different random samples
	if BOOTrep~=paramsOUT.BOOTSTRAP_Navgs+1
		%% Window spike times (ignore 1st 50 ms, and cut to shorter of 2 stim)
		[ST_A_plus,NumDrivenSpikes(1,1)]=windowSTs(SpikeTrains{1,1}(BOOTinds{1,1}{BOOTrep}), ...
			paramsOUT.SCC_onsetIGNORE_sec,stimdur_sec,paramsOUT.MAXspikes);
		[ST_A_minus,NumDrivenSpikes(1,2)]=windowSTs(SpikeTrains{1,2}(BOOTinds{1,2}{BOOTrep}), ...
			paramsOUT.SCC_onsetIGNORE_sec,stimdur_sec,paramsOUT.MAXspikes);
		[ST_B_plus,NumDrivenSpikes(2,1)]=windowSTs(SpikeTrains{2,1}(BOOTinds{2,1}{BOOTrep}), ...
			paramsOUT.SCC_onsetIGNORE_sec,stimdur_sec,paramsOUT.MAXspikes);
		[ST_B_minus,NumDrivenSpikes(2,2)]=windowSTs(SpikeTrains{2,2}(BOOTinds{2,2}{BOOTrep}), ...
			paramsOUT.SCC_onsetIGNORE_sec,stimdur_sec,paramsOUT.MAXspikes);
		paramsOUT.SCCdur_sec=stimdur_sec-paramsOUT.SCC_onsetIGNORE_sec;

		NumReps(1,1)=length(ST_A_plus); NumReps(1,2)=length(ST_A_minus); NumReps(2,1)=length(ST_B_plus); NumReps(2,2)=length(ST_B_minus);
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		disp(sprintf('NUMBER OF SPIKES (A_plus, A_minus, B_plus, B_minus) = %d, %d, %d, %d',NumDrivenSpikes(1,1),NumDrivenSpikes(1,2), ...
			NumDrivenSpikes(2,1),NumDrivenSpikes(2,2)))
		disp(sprintf('NUMBER OF REPS (A_plus, A_minus, B_plus, B_minus) = %d, %d, %d, %d',NumReps(1,1),NumReps(1,2),NumReps(2,1),NumReps(2,2)))

		disp(sprintf('   ... Computing SACs/SCCs for AN spikes: REP %d of %d ...',BOOTrep,paramsOUT.BOOTSTRAP_Navgs))
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

	else	% For LAST loop index, take AVG of all functions for computations
		% AVG everything! Then run rest of code ASIS!!
		% currently it has the last REP already
		for BOOTrep2=1:paramsOUT.BOOTSTRAP_Navgs-1
			SACSCCs.SAC_A_avg=SACSCCs.SAC_A_avg+SACSCCfunctions{BOOTrep2}.SAC_A_avg;
			SACSCCs.SAC_B_avg=SACSCCs.SAC_B_avg+SACSCCfunctions{BOOTrep2}.SAC_B_avg;
			SACSCCs.SCC_AB_avg=SACSCCs.SCC_AB_avg+SACSCCfunctions{BOOTrep2}.SCC_AB_avg;
			SACSCCs.XpAC_A_avg=SACSCCs.XpAC_A_avg+SACSCCfunctions{BOOTrep2}.XpAC_A_avg;
			SACSCCs.XpAC_B_avg=SACSCCs.XpAC_B_avg+SACSCCfunctions{BOOTrep2}.XpAC_B_avg;
			SACSCCs.XpCC_AB_avg=SACSCCs.XpCC_AB_avg+SACSCCfunctions{BOOTrep2}.XpCC_AB_avg;
			SACSCCs.SUMCOR_A=SACSCCs.SUMCOR_A+SACSCCfunctions{BOOTrep2}.SUMCOR_A;
			SACSCCs.SUMCOR_B=SACSCCs.SUMCOR_B+SACSCCfunctions{BOOTrep2}.SUMCOR_B;
			SACSCCs.SUMCOR_AB=SACSCCs.SUMCOR_AB+SACSCCfunctions{BOOTrep2}.SUMCOR_AB;
			SACSCCs.SUMCORadj_A=SACSCCs.SUMCORadj_A+SACSCCfunctions{BOOTrep2}.SUMCORadj_A;
			SACSCCs.SUMCORadj_B=SACSCCs.SUMCORadj_B+SACSCCfunctions{BOOTrep2}.SUMCORadj_B;
			SACSCCs.SUMCORadj_AB=SACSCCs.SUMCORadj_AB+SACSCCfunctions{BOOTrep2}.SUMCORadj_AB;
			SACSCCs.DIFCOR_A=SACSCCs.DIFCOR_A+SACSCCfunctions{BOOTrep2}.DIFCOR_A;
			SACSCCs.DIFCOR_B=SACSCCs.DIFCOR_B+SACSCCfunctions{BOOTrep2}.DIFCOR_B;
			SACSCCs.DIFCOR_AB=SACSCCs.DIFCOR_AB+SACSCCfunctions{BOOTrep2}.DIFCOR_AB;
			SACSCCs.PSDsc_A=SACSCCs.PSDsc_A+SACSCCfunctions{BOOTrep2}.PSDsc_A;
			SACSCCs.PSDsc_B=SACSCCs.PSDsc_B+SACSCCfunctions{BOOTrep2}.PSDsc_B;
			SACSCCs.CSDsc_AB=SACSCCs.CSDsc_AB+SACSCCfunctions{BOOTrep2}.CSDsc_AB;
			% RAND
			SACSCCs_rand.PSDsc_A=SACSCCs_rand.PSDsc_A+SACSCCfunctions{BOOTrep2}.rand.PSDsc_A;
			SACSCCs_rand.PSDsc_B=SACSCCs_rand.PSDsc_B+SACSCCfunctions{BOOTrep2}.rand.PSDsc_B;
			SACSCCs_rand.CSDsc_AB=SACSCCs_rand.CSDsc_AB+SACSCCfunctions{BOOTrep2}.rand.CSDsc_AB;
		end
		SACSCCs.SAC_A_avg=SACSCCs.SAC_A_avg/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs.SAC_B_avg=SACSCCs.SAC_B_avg/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs.SCC_AB_avg=SACSCCs.SCC_AB_avg/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs.XpAC_A_avg=SACSCCs.XpAC_A_avg/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs.XpAC_B_avg=SACSCCs.XpAC_B_avg/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs.XpCC_AB_avg=SACSCCs.XpCC_AB_avg/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs.SUMCOR_A=SACSCCs.SUMCOR_A/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs.SUMCOR_B=SACSCCs.SUMCOR_B/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs.SUMCOR_AB=SACSCCs.SUMCOR_AB/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs.SUMCORadj_A=SACSCCs.SUMCORadj_A/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs.SUMCORadj_B=SACSCCs.SUMCORadj_B/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs.SUMCORadj_AB=SACSCCs.SUMCORadj_AB/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs.DIFCOR_A=SACSCCs.DIFCOR_A/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs.DIFCOR_B=SACSCCs.DIFCOR_B/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs.DIFCOR_AB=SACSCCs.DIFCOR_AB/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs.PSDsc_A=SACSCCs.PSDsc_A/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs.PSDsc_B=SACSCCs.PSDsc_B/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs.CSDsc_AB=SACSCCs.CSDsc_AB/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs_rand.PSDsc_A=SACSCCs_rand.PSDsc_A/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs_rand.PSDsc_B=SACSCCs_rand.PSDsc_B/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs_rand.CSDsc_AB=SACSCCs_rand.CSDsc_AB/paramsOUT.BOOTSTRAP_Navgs;
		
		NumDrivenSpikes=zeros(2);
		AvgRate_sps=zeros(2);
		for BOOTrep3=1:paramsOUT.BOOTSTRAP_Navgs
			NumDrivenSpikes=NumDrivenSpikes+SACSCCmetrics{BOOTrep3}.NumDrivenSpikes/paramsOUT.BOOTSTRAP_Navgs;
			AvgRate_sps=AvgRate_sps+SACSCCmetrics{BOOTrep3}.AvgRate_sps/paramsOUT.BOOTSTRAP_Navgs;
		end
	end


	%%%%%%%%%%%%%%%%%%%%%
	%% PSD/CSD summations - Compute summed energy in ENVELOPE spectral
	%% densities over various frequency ranges
	%%%%%%%%%%%%%%%%%%%%%
	% create complete list of LH frequencies over which sums are needed for 
	% SCpeak_adjusted and/or CCCenv
	NYQ_Hz=0.5*(1/paramsOUT.DELAYbinwidth_sec);
	% default list
	default_PSD_LHfreqs_Hz=[0 NYQ_Hz; 0 CF_Hz; 0 300; 10 300; 10 CF_Hz; 0 50; 0 100; 0 64; 0 150];
	% add user-passed list, and create unique list
	if isfield(paramsOUT,'PSD_LHfreqs_Hz')
		paramsOUT.PSD_LHfreqs_Hz=unique([default_PSD_LHfreqs_Hz; paramsOUT.PSD_LHfreqs_Hz],'rows');
	else
		paramsOUT.PSD_LHfreqs_Hz=default_PSD_LHfreqs_Hz;
	end
	
	% find INDs in freqVEC_Hz for relevant cutoffs
	for i=1:size(paramsOUT.PSD_LHfreqs_Hz,1)
		[y,PSD_LHfreqs_inds(i,1)]=min(abs(SACSCCs.freqVEC-paramsOUT.PSD_LHfreqs_Hz(i,1)));
		[y,PSD_LHfreqs_inds(i,2)]=min(abs(SACSCCs.freqVEC-paramsOUT.PSD_LHfreqs_Hz(i,2)));
	end

	% Compute all sums for all spectral densities
	for i=1:size(paramsOUT.PSD_LHfreqs_Hz,1)
		SACSCCmetrics{BOOTrep}.sums.sumPSD_A(i) = sum(SACSCCs.PSDsc_A(PSD_LHfreqs_inds(i,1):PSD_LHfreqs_inds(i,2)));
		if ~isnan(NumDrivenSpikes(2,1))
			SACSCCmetrics{BOOTrep}.sums.sumPSD_B(i) = sum(SACSCCs.PSDsc_B(PSD_LHfreqs_inds(i,1):PSD_LHfreqs_inds(i,2)));
			SACSCCmetrics{BOOTrep}.sums.sumCSD_AB(i) = sum(SACSCCs.CSDsc_AB(PSD_LHfreqs_inds(i,1):PSD_LHfreqs_inds(i,2)));
			SACSCCmetrics{BOOTrep}.sums.sumPSDrand_A(i) = sum(SACSCCs_rand.PSDsc_A(PSD_LHfreqs_inds(i,1):PSD_LHfreqs_inds(i,2)));
			SACSCCmetrics{BOOTrep}.sums.sumPSDrand_B(i) = sum(SACSCCs_rand.PSDsc_B(PSD_LHfreqs_inds(i,1):PSD_LHfreqs_inds(i,2)));
			SACSCCmetrics{BOOTrep}.sums.sumCSDrand_AB(i) = sum(SACSCCs_rand.CSDsc_AB(PSD_LHfreqs_inds(i,1):PSD_LHfreqs_inds(i,2)));
		end
	end
	SACSCCmetrics{BOOTrep}.sums.PSD_LHfreqs_Hz=paramsOUT.PSD_LHfreqs_Hz;
		
	%%%%%%%%%%%%%%%%%%%%
	%% Compute metrics
	%%%%%%%%%%%%%%%%%%%%
	%% SAC/DC/SC Peak Heights
	% SAC peaks - this is CI from Joris et al 2006 (HR) [don't subtract 1]
	SACSCCmetrics{BOOTrep}.SACpeak_A=max(SACSCCs.SAC_A_avg);
	if ~isnan(NumDrivenSpikes(2,1))
		SACSCCmetrics{BOOTrep}.SACpeak_B=max(SACSCCs.SAC_B_avg);
		SACSCCmetrics{BOOTrep}.SCCpeak_AB=max(SACSCCs.SCC_AB_avg);
	end
	% DIFCOR peaks
	SACSCCmetrics{BOOTrep}.DCpeak_A=max(SACSCCs.DIFCOR_A);
	if ~isnan(NumDrivenSpikes(2,1))
		SACSCCmetrics{BOOTrep}.DCpeak_B=max(SACSCCs.DIFCOR_B);
		SACSCCmetrics{BOOTrep}.DCpeak_AB=max(SACSCCs.DIFCOR_AB);
	end
	% SUMCOR peaks (don't subtract 1; Louage et al 2004)
	% compute three default SCpeaks
	% 1) raw peaks
	SACSCCmetrics{BOOTrep}.SCpeaks_legend{1}='raw';
	SACSCCmetrics{BOOTrep}.SCpeaks_A(1)=max(SACSCCs.SUMCOR_A);
	if ~isnan(NumDrivenSpikes(2,1))
		SACSCCmetrics{BOOTrep}.SCpeaks_B(1)=max(SACSCCs.SUMCOR_B);
		SACSCCmetrics{BOOTrep}.SCpeaks_AB(1)=max(SACSCCs.SUMCOR_AB);
	end
	% 2) Adjusted SCpeak (Eq. 2): Adj SCpeak = rawSCpeak*sum[0,CF]/sum[0,Fs/2]
	%% USE THIS ONE: 
	SACSCCmetrics{BOOTrep}.SCpeaks_legend{2}='adj: 0-CF';
	CFsum_index=find((paramsOUT.PSD_LHfreqs_Hz(:,1)==0)&(paramsOUT.PSD_LHfreqs_Hz(:,2)==CF_Hz));
	NYQsum_index=find((paramsOUT.PSD_LHfreqs_Hz(:,1)==0)&(paramsOUT.PSD_LHfreqs_Hz(:,2)==NYQ_Hz));
	SACSCCmetrics{BOOTrep}.SCpeaks_A(2)=1+(max(SACSCCs.SUMCOR_A)-1)* ...
		SACSCCmetrics{BOOTrep}.sums.sumPSD_A(CFsum_index)/SACSCCmetrics{BOOTrep}.sums.sumPSD_A(NYQsum_index);
	if ~isnan(NumDrivenSpikes(2,1))
		SACSCCmetrics{BOOTrep}.SCpeaks_B(2)=1+(max(SACSCCs.SUMCOR_B)-1)* ...
			SACSCCmetrics{BOOTrep}.sums.sumPSD_B(CFsum_index)/SACSCCmetrics{BOOTrep}.sums.sumPSD_B(NYQsum_index);
		SACSCCmetrics{BOOTrep}.SCpeaks_AB(2)=1+(max(SACSCCs.SUMCOR_AB)-1)* ...
			SACSCCmetrics{BOOTrep}.sums.sumCSD_AB(CFsum_index)/SACSCCmetrics{BOOTrep}.sums.sumCSD_AB(NYQsum_index);
	end
	% 3) alternative adjusted SC peak height [0,300]/[0,Fs/2]
	% just for FIG 2 in NM paper
	SACSCCmetrics{BOOTrep}.SCpeaks_legend{3}='adj: 0-300';
	ENVsum_index=find((paramsOUT.PSD_LHfreqs_Hz(:,1)==0)&(paramsOUT.PSD_LHfreqs_Hz(:,2)==300));
	SACSCCmetrics{BOOTrep}.SCpeaks_A(3)=1+(max(SACSCCs.SUMCOR_A)-1)* ...
		SACSCCmetrics{BOOTrep}.sums.sumPSD_A(ENVsum_index)/SACSCCmetrics{BOOTrep}.sums.sumPSD_A(NYQsum_index);
	if ~isnan(NumDrivenSpikes(2,1))
		SACSCCmetrics{BOOTrep}.SCpeaks_B(3)=1+(max(SACSCCs.SUMCOR_B)-1)* ...
			SACSCCmetrics{BOOTrep}.sums.sumPSD_B(ENVsum_index)/SACSCCmetrics{BOOTrep}.sums.sumPSD_B(NYQsum_index);
		SACSCCmetrics{BOOTrep}.SCpeaks_AB(3)=1+(max(SACSCCs.SUMCOR_AB)-1)* ...
			SACSCCmetrics{BOOTrep}.sums.sumCSD_AB(ENVsum_index)/SACSCCmetrics{BOOTrep}.sums.sumCSD_AB(NYQsum_index);
	end
	% 4) raw peaks
	SACSCCmetrics{BOOTrep}.SCpeaks_legend{4}='IFFTraw';
	SACSCCmetrics{BOOTrep}.SCpeaks_A(4)=max(SACSCCs.SUMCORadj_A);
	if ~isnan(NumDrivenSpikes(2,1))
		SACSCCmetrics{BOOTrep}.SCpeaks_B(4)=max(SACSCCs.SUMCORadj_B);
		SACSCCmetrics{BOOTrep}.SCpeaks_AB(4)=max(SACSCCs.SUMCORadj_AB);
	end
	
	%% Neural Cross Correlation Coefficients
	if ~isnan(NumDrivenSpikes(2,1))
		%% CCCtfs
		if (SACSCCmetrics{BOOTrep}.DCpeak_A>=paramsOUT.minDCpeak_CCCtfs)&(SACSCCmetrics{BOOTrep}.DCpeak_B>=paramsOUT.minDCpeak_CCCtfs)
			SACSCCmetrics{BOOTrep}.CCCtfs = SACSCCmetrics{BOOTrep}.DCpeak_AB/(sqrt(SACSCCmetrics{BOOTrep}.DCpeak_A*SACSCCmetrics{BOOTrep}.DCpeak_B));
		else
			SACSCCmetrics{BOOTrep}.CCCtfs = NaN;
		end
		%% CCCenv
		% first, go through whole LHfreq list, for subBIAS (subtrat bias) and
		% withBIAS (Don't subtract bias)
		for i=1:size(paramsOUT.PSD_LHfreqs_Hz,1)
			if paramsOUT.PSD_LHfreqs_Hz(i,2)==CF_Hz
				SACSCCmetrics{BOOTrep}.CCCenvs_legend{(i-1)*2+1}=sprintf('%.f-CF, subBIAS',paramsOUT.PSD_LHfreqs_Hz(i,1));
				SACSCCmetrics{BOOTrep}.CCCenvs_legend{(i-1)*2+2}=sprintf('%.f-CF, withBIAS',paramsOUT.PSD_LHfreqs_Hz(i,1));
			elseif paramsOUT.PSD_LHfreqs_Hz(i,2)==NYQ_Hz
				SACSCCmetrics{BOOTrep}.CCCenvs_legend{(i-1)*2+1}=sprintf('%.f-NYQ, subBIAS',paramsOUT.PSD_LHfreqs_Hz(i,1));
				SACSCCmetrics{BOOTrep}.CCCenvs_legend{(i-1)*2+2}=sprintf('%.f-NYQ, withBIAS',paramsOUT.PSD_LHfreqs_Hz(i,1));
			else
				SACSCCmetrics{BOOTrep}.CCCenvs_legend{(i-1)*2+1}=sprintf('%.f-%.f, subBIAS',paramsOUT.PSD_LHfreqs_Hz(i,1),paramsOUT.PSD_LHfreqs_Hz(i,2));
				SACSCCmetrics{BOOTrep}.CCCenvs_legend{(i-1)*2+2}=sprintf('%.f-%.f, withBIAS',paramsOUT.PSD_LHfreqs_Hz(i,1),paramsOUT.PSD_LHfreqs_Hz(i,2));
			end
			% subBIAS
			SACSCCmetrics{BOOTrep}.CCCenvs((i-1)*2+1) = max([0 (SACSCCmetrics{BOOTrep}.sums.sumCSD_AB(i)-SACSCCmetrics{BOOTrep}.sums.sumCSDrand_AB(i)) ])/ ...
				sqrt(max([0 (SACSCCmetrics{BOOTrep}.sums.sumPSD_A(i)-SACSCCmetrics{BOOTrep}.sums.sumPSDrand_A(i)) ])* ...
				max([0 (SACSCCmetrics{BOOTrep}.sums.sumPSD_B(i)-SACSCCmetrics{BOOTrep}.sums.sumPSDrand_B(i)) ]));
			% withBIAS
			SACSCCmetrics{BOOTrep}.CCCenvs((i-1)*2+2) = SACSCCmetrics{BOOTrep}.sums.sumCSD_AB(i)/ ...
				sqrt(SACSCCmetrics{BOOTrep}.sums.sumPSD_A(i)*SACSCCmetrics{BOOTrep}.sums.sumPSD_B(i));
		end
		% using raw SCpeaks
		RAWSCindex=find(strcmp(SACSCCmetrics{BOOTrep}.SCpeaks_legend,'raw'));
		SACSCCmetrics{BOOTrep}.CCCenvs_legend{end+1}='rawSC';
		SACSCCmetrics{BOOTrep}.CCCenvs(end+1) = (SACSCCmetrics{BOOTrep}.SCpeaks_AB(RAWSCindex)-1)/ ...
			(sqrt((SACSCCmetrics{BOOTrep}.SCpeaks_A(RAWSCindex)-1)*(SACSCCmetrics{BOOTrep}.SCpeaks_B(RAWSCindex)-1)));
		% using adj SCpeaks
		ADJSCindex=find(strcmp(SACSCCmetrics{BOOTrep}.SCpeaks_legend,'adj: 0-CF'));
		SACSCCmetrics{BOOTrep}.CCCenvs_legend{end+1}='adjSC';
		SACSCCmetrics{BOOTrep}.CCCenvs(end+1) = (SACSCCmetrics{BOOTrep}.SCpeaks_AB(ADJSCindex)-1)/ ...
			(sqrt((SACSCCmetrics{BOOTrep}.SCpeaks_A(ADJSCindex)-1)*(SACSCCmetrics{BOOTrep}.SCpeaks_B(ADJSCindex)-1)));
		% using IFFTraw SCpeaks
		IFFTSCindex=find(strcmp(SACSCCmetrics{BOOTrep}.SCpeaks_legend,'IFFTraw'));
		SACSCCmetrics{BOOTrep}.CCCenvs_legend{end+1}='IFFTrawSC';
		SACSCCmetrics{BOOTrep}.CCCenvs(end+1) = (SACSCCmetrics{BOOTrep}.SCpeaks_AB(IFFTSCindex)-1)/ ...
			(sqrt((SACSCCmetrics{BOOTrep}.SCpeaks_A(IFFTSCindex)-1)*(SACSCCmetrics{BOOTrep}.SCpeaks_B(IFFTSCindex)-1)));
		
		%% Characteristic Delays
		%% NEED TO IMPROVE - find largest, unless there is a "close-to-largest" that
		% is closer to zero!
		% probably use some criterion for 2nd largest peak (if with 5%)???
		SACSCCmetrics{BOOTrep}.CDscc_usec=findCD_SCC(SACSCCs.SCC_AB_avg,SACSCCs.delays_usec);
		SACSCCmetrics{BOOTrep}.CDenv_usec=findCD_SCC(SACSCCs.SUMCORadj_AB,SACSCCs.delays_usec);  % use IFFTadjusted SUMCOR!
		SACSCCmetrics{BOOTrep}.CDtfs_usec=findCD_SCC(SACSCCs.DIFCOR_AB,SACSCCs.delays_usec);
	end

	%% Book-keeping
	SACSCCfunctions{BOOTrep}=SACSCCs;
	if ~isnan(NumDrivenSpikes(2,1))
		SACSCCfunctions{BOOTrep}.rand.PSDsc_A=SACSCCs_rand.PSDsc_A;
		SACSCCfunctions{BOOTrep}.rand.PSDsc_B=SACSCCs_rand.PSDsc_B;
		SACSCCfunctions{BOOTrep}.rand.CSDsc_AB=SACSCCs_rand.CSDsc_AB;
	end

	SACSCCmetrics{BOOTrep}.NumDrivenSpikes=NumDrivenSpikes;
	SACSCCmetrics{BOOTrep}.AvgRate_sps=AvgRate_sps;

% 	if PLOT_15panel
% 		plot_CCCanal_0(SACSCCfunctions{BOOTrep},SACSCCmetrics{BOOTrep},paramsOUT)
% 	end

end


%% %%%%%%%%%%%%%%% TODO 8/5/08
% *1) finish this for TOUSE SCpeaks CCCenv  BOOKKEEPING ...
%    - go ahead and compute VEC and AVG/STD for ALL CCCenv computations,
%    then just access LATER!!!!
% *2) fix plot_CCCanal_3 as in plot_0
% *3) test plots with user defined TOUSEs after computations
% *- setup general to output various options
% *4) generalize GEN 2 conditions
%
%% 5) GANESH PROOF: verify vs CCCanal_3a
%% 7) test with real data - BBN ABA, CHIM 1, 16 
%   - decide how to do comps!
%      - bootstrap
%      - freq range (BBN, CHIM)
%      - number of spikes 
% 6) test Fig. 2
%% 7) GO OVER WITH GANESH
% - him proof it all
%% - setup BML code

%% TALLY AVG
CCCtfsVEC=zeros(1,paramsOUT.BOOTSTRAP_Navgs);
for BOOTrep4=1:paramsOUT.BOOTSTRAP_Navgs
	CCCtfsVEC(BOOTrep4)=SACSCCmetrics{BOOTrep4}.CCCtfs;
end
CCCtfsAVG=mean(CCCtfsVEC);
CCCtfsSTD=std(CCCtfsVEC);
CCCenvVECs=cell(size(SACSCCmetrics{1}.CCCenvs));
CCCenvAVGs=zeros(size(SACSCCmetrics{1}.CCCenvs));
CCCenvSTDs=zeros(size(SACSCCmetrics{1}.CCCenvs));
for i=1:length(CCCenvVECs)
	CCCenvVECs{i}=zeros(1,paramsOUT.BOOTSTRAP_Navgs);
	for BOOTrep4=1:paramsOUT.BOOTSTRAP_Navgs
		CCCenvVECs{i}(BOOTrep4)=SACSCCmetrics{BOOTrep4}.CCCenvs(i);
	end
	CCCenvAVGs(i)=mean(CCCenvVECs{i});
	CCCenvSTDs(i)=std(CCCenvVECs{i});
end
SACSCCmetrics{paramsOUT.BOOTSTRAP_Navgs+1}.CCCenvVECs=CCCenvVECs;
SACSCCmetrics{paramsOUT.BOOTSTRAP_Navgs+1}.CCCenvAVGs=CCCenvAVGs;
SACSCCmetrics{paramsOUT.BOOTSTRAP_Navgs+1}.CCCenvSTDs=CCCenvSTDs;
SACSCCmetrics{paramsOUT.BOOTSTRAP_Navgs+1}.CCCtfsVEC=CCCtfsVEC;
SACSCCmetrics{paramsOUT.BOOTSTRAP_Navgs+1}.CCCtfsAVG=CCCtfsAVG;
SACSCCmetrics{paramsOUT.BOOTSTRAP_Navgs+1}.CCCtfsSTD=CCCtfsSTD;

if PLOT_15panel
	plot_CCCanal_3(SACSCCfunctions,SACSCCmetrics,paramsOUT)
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

%% Remove TFS artifact (centered at 2*CF) from SUMCOR
% zero out above CF
FFTtemp=fft((SUMCOR_A-1),paramsOUT.Nfft_psd);
freqVEC=(0:length(FFTtemp)-1)/length(FFTtemp)*paramsOUT.Nfft_psd;
[y,CF_index]=min(abs(freqVEC-paramsOUT.CF_A_Hz)); % use CF_A
FFTadj=zeros(size(FFTtemp));
FFTadj(1:CF_index)=FFTtemp(1:CF_index);
FFTadj((length(FFTtemp)-CF_index+1):end)=FFTtemp((length(FFTtemp)-CF_index+1):end); %keep negative freqs
adjSC=ifft(FFTadj)+1;
SUMCORadj_A=real(adjSC(1:length(SUMCOR_A)));

%% Compute ENVELOPE POWER SPECTRAL DENSITY for Condition A
PSDsc_A=abs(FFTtemp);

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

	%% Remove TFS artifact (centered at 2*CF) from SUMCOR
	% zero out above CF
	FFTtemp=fft((SUMCOR_B-1),paramsOUT.Nfft_psd);
	[y,CF_index]=min(abs(freqVEC-paramsOUT.CF_B_Hz)); % use CF_B
	FFTadj=zeros(size(FFTtemp));
	FFTadj(1:CF_index)=FFTtemp(1:CF_index);
	FFTadj((length(FFTtemp)-CF_index+1):end)=FFTtemp((length(FFTtemp)-CF_index+1):end); %keep negative freqs
	adjSC=ifft(FFTadj)+1;
	SUMCORadj_B=real(adjSC(1:length(SUMCOR_B)));

	%% Compute ENVELOPE POWER SPECTRAL DENSITY for Condition B
	PSDsc_B=abs(FFTtemp);
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

	%% Remove TFS artifact (centered at 2*CF) from SUMCOR
	% zero out above CF
	FFTtemp=fft((SUMCOR_AB-1),paramsOUT.Nfft_psd);
	[y,CF_index]=min(abs(freqVEC-paramsOUT.SACSCC_CF_Hz)); % use min([CF_A CF_B])
	FFTadj=zeros(size(FFTtemp));
	FFTadj(1:CF_index)=FFTtemp(1:CF_index);
	FFTadj((length(FFTtemp)-CF_index+1):end)=FFTtemp((length(FFTtemp)-CF_index+1):end); %keep negative freqs
	adjSC=ifft(FFTadj)+1;
	SUMCORadj_AB=real(adjSC(1:length(SUMCOR_AB)));

	%% Compute ENVELOPE CROSS SPECTRAL DENSITY for Condition AB
	CSDsc_AB=abs(FFTtemp);
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
	SUMCORadj_B=[];
	PSDsc_B=[];
	SCC_AB_plus=[];
	SCC_AB_minus=[];
	SCC_AB_avg=[];
	XpCC_AB_plus=[];
	XpCC_AB_minus=[];
	XpCC_AB_avg=[];
	DIFCOR_AB=[];
	SUMCOR_AB=[];
	SUMCORadj_AB=[];
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
	'SUMCOR_B',SUMCOR_B,'SUMCOR_AB',SUMCOR_AB,'SUMCORadj_A',SUMCORadj_A,'SUMCORadj_B',SUMCORadj_B,'SUMCORadj_AB',SUMCORadj_AB, ...
	'DIFCOR_A',DIFCOR_A,'DIFCOR_B',DIFCOR_B,'DIFCOR_AB',DIFCOR_AB,'PSDsc_A',PSDsc_A,'PSDsc_B',PSDsc_B,'CSDsc_AB',CSDsc_AB,'freqVEC',freqVEC);

return;

