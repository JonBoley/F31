function PATCH_100408_CCCanal_3(ORIGfilename)
% File: PATCH_100408_CCCanal_3.m
% M. Heinz
% Oct 4, 2008 - re-compute CD for Delta_CF with find_CDpos_SCC - to mak
% esure CD is > 0 (USED only for 1000Hz CF Delta-CD vs LEVEL
%
%
% Sept 26, 2008
%
% Modified: from CCCanal_3, to PATCH data already run with old CCCanal_3
% modifications are:
% 1) PSDs/CSD for AN spikes for AVG_BOOTSTRAP condition are no longer taken
% as AVG of 5 PSDs, but rather the PSD of the AVG SCs.  Thought this might
% reduce overestimate of CCCenv in AA conditions (which have been above 1,
% especially for the AN data
% 2) IFFT-based SC computations are expanded to compute SCpeaks at several
% specified delays, in addition to just taking the max: Delays = 0, CD, and
% user-specified delays.  Note: CCCs based on IFFT_SC peaks use delay of 0
% for ACFs, and only specified delay for CCF (AB) condition.
%
% Basics of PATCH: 
%  -- only changes to SACSCCfunctions: recompute AVG PSDs/CSD
%  -- completely recommpute SACSCCmetrics!
%  -- save as NEWFileName (w, '_PATCH100408' behind ORIGFileName)

if ~exist('ORIGfilename','var'),    ORIGfilename='BBNAB_CF422_SCCs';   end
NEWfilename=[ORIGfilename '_PATCH100408'];

%% LOAD existing CCCanal Data
disp(sprintf('LOADING EXISTING CCCanal_3 file: %s.mat',ORIGfilename))
eval(['load ' ORIGfilename '.mat'])

plot_CCCanal_3p(SACSCCfunctions,SACSCCmetrics,paramsOUT)
set(gcf,'Name',sprintf('Pre-PATCH: %s',ORIGfilename))

%% Clear SACSCCmetrics (save old for archival purposes)
SACSCCmetrics_prePATCH100408=SACSCCmetrics;
clear SACSCCmetrics
paramsOUT_prePATCH100408=paramsOUT;

%% Update paramsOUT %%%%%%%
% Default SACpeak, DCpeak, SCpeak and CCCenvs to use
paramsOUT.SACpeak_TOUSE='SACpeak_0';
paramsOUT.DCpeak_TOUSE='DCpeak_0';
paramsOUT.CCCtfs_TOUSE='DCpeak_0';
paramsOUT.SCpeak_TOUSE='IFFTraw_0';
paramsOUT.CCCenv_TOUSE='IFFTrawSC_0';


%%%%%%%%%%%%%%%%%%%%%%%
%% Recompute PSDs/CSD based on AVG SCs
%%%%%%%%%%%%%%%%%%%%%%%
SACSCCfunctions{end}.PSDsc_A=abs(fft((SACSCCfunctions{end}.SUMCOR_A-1),paramsOUT.Nfft_psd));
SACSCCfunctions{end}.PSDsc_B=abs(fft((SACSCCfunctions{end}.SUMCOR_B-1),paramsOUT.Nfft_psd));
SACSCCfunctions{end}.CSDsc_AB=abs(fft((SACSCCfunctions{end}.SUMCOR_AB-1),paramsOUT.Nfft_psd));

% plot_CCCanal_3p(SACSCCfunctions,SACSCCmetrics_prePATCH100408,paramsOUT_prePATCH100408)
% set(gcf,'Name',sprintf('PATCH_PSD: %s',ORIGfilename))


%% Recompute SACSCCmetrics for all BOOTreps
for BOOTrep=1:paramsOUT.BOOTSTRAP_Navgs+1

	% just to keep code the same as from CCCanal_3
	SACSCCs=SACSCCfunctions{BOOTrep};	
	SACSCCs_rand=SACSCCfunctions{BOOTrep}.rand;
	NumDrivenSpikes=SACSCCmetrics_prePATCH100408{BOOTrep}.NumDrivenSpikes;
	CF_Hz=paramsOUT.SACSCC_CF_Hz;
	
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

	%% Characteristic Delays
	%% NEED TO IMPROVE - find largest, unless there is a "close-to-largest" that
	%% is closer to zero!
	% probably use some criterion for 2nd largest peak (if with 5%)???
	%
	% CD > 0 enforced HERE
	%
	SACSCCmetrics{BOOTrep}.CDscc_usec=findCDpos_SCC(SACSCCs.SCC_AB_avg,SACSCCs.delays_usec);
	SACSCCmetrics{BOOTrep}.CDenv_usec=findCDpos_SCC(SACSCCs.SUMCORadj_AB,SACSCCs.delays_usec);  % use IFFTadjusted SUMCOR!
	SACSCCmetrics{BOOTrep}.CDtfs_usec=findCDpos_SCC(SACSCCs.DIFCOR_AB,SACSCCs.delays_usec);

	
	%% SAC/DC/SC Peak Heights
	%%%%%%%%%
	%% SAC peaks - this is CI from Joris et al 2006 (HR) [don't subtract 1]
	% 1) SACpeak (pure MAX)
	SACSCCmetrics{BOOTrep}.SACpeaks_legend{1}='SACpeak_max';
	SACSCCmetrics{BOOTrep}.SACpeak_A(1)=max(SACSCCs.SAC_A_avg);
	if ~isnan(NumDrivenSpikes(2,1))
		SACSCCmetrics{BOOTrep}.SACpeak_B(1)=max(SACSCCs.SAC_B_avg);
		SACSCCmetrics{BOOTrep}.SCCpeak_AB(1)=max(SACSCCs.SCC_AB_avg);
	end
	% 2) SACpeak (0 delay)
	SACSCCmetrics{BOOTrep}.SACpeaks_legend{2}='SACpeak_0';
	INDEX_0=find(SACSCCs.delays_usec==0);
	SACSCCmetrics{BOOTrep}.SACpeak_A(2)=SACSCCs.SAC_A_avg(INDEX_0);
	if ~isnan(NumDrivenSpikes(2,1))
		SACSCCmetrics{BOOTrep}.SACpeak_B(2)=SACSCCs.SAC_B_avg(INDEX_0);
		SACSCCmetrics{BOOTrep}.SCCpeak_AB(2)=SACSCCs.SCC_AB_avg(INDEX_0);
	end
	% 3) SACpeak (CD delay)
	SACSCCmetrics{BOOTrep}.SACpeaks_legend{3}='SACpeak_CD';
	INDEX_CD=find(SACSCCs.delays_usec==SACSCCmetrics{BOOTrep}.CDscc_usec);
	SACSCCmetrics{BOOTrep}.SACpeak_A(3)=SACSCCs.SAC_A_avg(INDEX_CD);
	if ~isnan(NumDrivenSpikes(2,1))
		SACSCCmetrics{BOOTrep}.SACpeak_B(3)=SACSCCs.SAC_B_avg(INDEX_CD);
		SACSCCmetrics{BOOTrep}.SCCpeak_AB(3)=SACSCCs.SCC_AB_avg(INDEX_CD);
	end
	%% User-passed list of delays to compute
	if isfield(paramsOUT,'UserDelays_usec')
		for i=1:length(paramsOUT.UserDelays_usec)
			% 4+) SACpeak (user delays)
			SACSCCmetrics{BOOTrep}.SACpeaks_legend{end+1}=sprintf('SACpeak_%d',round(paramsOUT.UserDelays_usec(i)));
			[yTEMP,INDEX_user]=min(abs(SACSCCs.delays_usec-paramsOUT.UserDelays_usec(i))); % find closest delay
			SACSCCmetrics{BOOTrep}.SACpeak_A(end+1)=SACSCCs.SAC_A_avg(INDEX_user);
			if ~isnan(NumDrivenSpikes(2,1))
				SACSCCmetrics{BOOTrep}.SACpeak_B(end+1)=SACSCCs.SAC_B_avg(INDEX_user);
				SACSCCmetrics{BOOTrep}.SCCpeak_AB(end+1)=SACSCCs.SCC_AB_avg(INDEX_user);
			end
		end
	end

	%%%%%%%%%
	%% DC peaks - 
	% 1) DCpeak (pure MAX)
	SACSCCmetrics{BOOTrep}.DCpeaks_legend{1}='DCpeak_max';
	SACSCCmetrics{BOOTrep}.DCpeak_A(1)=max(SACSCCs.DIFCOR_A);
	if ~isnan(NumDrivenSpikes(2,1))
		SACSCCmetrics{BOOTrep}.DCpeak_B(1)=max(SACSCCs.DIFCOR_B);
		SACSCCmetrics{BOOTrep}.DCpeak_AB(1)=max(SACSCCs.DIFCOR_AB);
	end
	% 2) DCpeak (0 delay)
	SACSCCmetrics{BOOTrep}.DCpeaks_legend{2}='DCpeak_0';
	SACSCCmetrics{BOOTrep}.DCpeak_A(2)=SACSCCs.DIFCOR_A(INDEX_0);
	if ~isnan(NumDrivenSpikes(2,1))
		SACSCCmetrics{BOOTrep}.DCpeak_B(2)=SACSCCs.DIFCOR_B(INDEX_0);
		SACSCCmetrics{BOOTrep}.DCpeak_AB(2)=SACSCCs.DIFCOR_AB(INDEX_0);
	end
	% 3) DCpeak (CD delay)
	SACSCCmetrics{BOOTrep}.DCpeaks_legend{3}='DCpeak_CD';
	INDEX_CD=find(SACSCCs.delays_usec==SACSCCmetrics{BOOTrep}.CDtfs_usec);
	SACSCCmetrics{BOOTrep}.DCpeak_A(3)=SACSCCs.DIFCOR_A(INDEX_CD);
	if ~isnan(NumDrivenSpikes(2,1))
		SACSCCmetrics{BOOTrep}.DCpeak_B(3)=SACSCCs.DIFCOR_B(INDEX_CD);
		SACSCCmetrics{BOOTrep}.DCpeak_AB(3)=SACSCCs.DIFCOR_AB(INDEX_CD);
	end
	%% User-passed list of delays to compute
	if isfield(paramsOUT,'UserDelays_usec')
		for i=1:length(paramsOUT.UserDelays_usec)
			% 4+) DCpeak (user delays)
			SACSCCmetrics{BOOTrep}.DCpeaks_legend{end+1}=sprintf('DCpeak_%d',round(paramsOUT.UserDelays_usec(i)));
			[yTEMP,INDEX_user]=min(abs(SACSCCs.delays_usec-paramsOUT.UserDelays_usec(i))); % find closest delay
			SACSCCmetrics{BOOTrep}.DCpeak_A(end+1)=SACSCCs.DIFCOR_A(INDEX_user);
			if ~isnan(NumDrivenSpikes(2,1))
				SACSCCmetrics{BOOTrep}.DCpeak_B(end+1)=SACSCCs.DIFCOR_B(INDEX_user);
				SACSCCmetrics{BOOTrep}.DCpeak_AB(end+1)=SACSCCs.DIFCOR_AB(INDEX_user);
			end
		end
	end

	%%%%%%%%%
	%% SUMCOR peaks (don't subtract 1; Louage et al 2004)
	% 1) raw peaks
	SACSCCmetrics{BOOTrep}.SCpeaks_legend{1}='raw';
	SACSCCmetrics{BOOTrep}.SCpeaks_A(1)=max(SACSCCs.SUMCOR_A);
	if ~isnan(NumDrivenSpikes(2,1))
		SACSCCmetrics{BOOTrep}.SCpeaks_B(1)=max(SACSCCs.SUMCOR_B);
		SACSCCmetrics{BOOTrep}.SCpeaks_AB(1)=max(SACSCCs.SUMCOR_AB);
	end
	% 2) Adjusted SCpeak (Eq. 2): Adj SCpeak = rawSCpeak*sum[0,CF]/sum[0,Fs/2]
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
	% 4) raw peaks of IFFTraw
	SACSCCmetrics{BOOTrep}.SCpeaks_legend{4}='IFFTraw';
	SACSCCmetrics{BOOTrep}.SCpeaks_A(4)=max(SACSCCs.SUMCORadj_A);
	if ~isnan(NumDrivenSpikes(2,1))
		SACSCCmetrics{BOOTrep}.SCpeaks_B(4)=max(SACSCCs.SUMCORadj_B);
		SACSCCmetrics{BOOTrep}.SCpeaks_AB(4)=max(SACSCCs.SUMCORadj_AB);
	end
	% 5) IFFTraw (0 delay)
	SACSCCmetrics{BOOTrep}.SCpeaks_legend{5}='IFFTraw_0';
	INDEX_0=find(SACSCCs.delays_usec==0);
	SACSCCmetrics{BOOTrep}.SCpeaks_A(5)=SACSCCs.SUMCORadj_A(INDEX_0);
	if ~isnan(NumDrivenSpikes(2,1))
		SACSCCmetrics{BOOTrep}.SCpeaks_B(5)=SACSCCs.SUMCORadj_B(INDEX_0);
		SACSCCmetrics{BOOTrep}.SCpeaks_AB(5)=SACSCCs.SUMCORadj_AB(INDEX_0);
	end
	% 6) IFFTraw (CD delay)
	SACSCCmetrics{BOOTrep}.SCpeaks_legend{end+1}='IFFTraw_CD';
	INDEX_CD=find(SACSCCs.delays_usec==SACSCCmetrics{BOOTrep}.CDenv_usec);
	SACSCCmetrics{BOOTrep}.SCpeaks_A(end+1)=SACSCCs.SUMCORadj_A(INDEX_CD);
	if ~isnan(NumDrivenSpikes(2,1))
		SACSCCmetrics{BOOTrep}.SCpeaks_B(end+1)=SACSCCs.SUMCORadj_B(INDEX_CD);
		SACSCCmetrics{BOOTrep}.SCpeaks_AB(end+1)=SACSCCs.SUMCORadj_AB(INDEX_CD);
	end
	%% User-passed list of delays to compute
	if isfield(paramsOUT,'UserDelays_usec')
		for i=1:length(paramsOUT.UserDelays_usec)
			% 7+) IFFTraw (user delays)
			SACSCCmetrics{BOOTrep}.SCpeaks_legend{end+1}=sprintf('IFFTraw_%d',round(paramsOUT.UserDelays_usec(i)));
			[yTEMP,INDEX_user]=min(abs(SACSCCs.delays_usec-paramsOUT.UserDelays_usec(i))); % find closest delay
			SACSCCmetrics{BOOTrep}.SCpeaks_A(end+1)=SACSCCs.SUMCORadj_A(INDEX_user);
			if ~isnan(NumDrivenSpikes(2,1))
				SACSCCmetrics{BOOTrep}.SCpeaks_B(end+1)=SACSCCs.SUMCORadj_B(INDEX_user);
				SACSCCmetrics{BOOTrep}.SCpeaks_AB(end+1)=SACSCCs.SUMCORadj_AB(INDEX_user);
			end
		end
	end

	%%%%%%%%%%%%%%%%%
	%% Neural Cross Correlation Coefficients
	%%%%%%%%%%%%%%%%%
	if ~isnan(NumDrivenSpikes(2,1))
		%% CCCtfs
		DC_max_index=find(strcmp(SACSCCmetrics{BOOTrep}.DCpeaks_legend,'DCpeak_max'));
		if (SACSCCmetrics{BOOTrep}.DCpeak_A(DC_max_index)>=paramsOUT.minDCpeak_CCCtfs) ...
				& (SACSCCmetrics{BOOTrep}.DCpeak_B(DC_max_index)>=paramsOUT.minDCpeak_CCCtfs)
			% 1) using DCpeak_max
			SACSCCmetrics{BOOTrep}.CCCtfs_legend{1}='DCpeak_max';
			SACSCCmetrics{BOOTrep}.CCCtfs(1) = SACSCCmetrics{BOOTrep}.DCpeak_AB(DC_max_index) ...
				/ (sqrt(SACSCCmetrics{BOOTrep}.DCpeak_A(DC_max_index)*SACSCCmetrics{BOOTrep}.DCpeak_B(DC_max_index)));
			% 2) using DCpeak_0
			DC_0_index=find(strcmp(SACSCCmetrics{BOOTrep}.DCpeaks_legend,'DCpeak_0'));
			SACSCCmetrics{BOOTrep}.CCCtfs_legend{2}='DCpeak_0';
			SACSCCmetrics{BOOTrep}.CCCtfs(2) = SACSCCmetrics{BOOTrep}.DCpeak_AB(DC_0_index) ...
				/ (sqrt(SACSCCmetrics{BOOTrep}.DCpeak_A(DC_0_index)*SACSCCmetrics{BOOTrep}.DCpeak_B(DC_0_index)));
			% 3) using DCpeak_CD
			DC_CD_index=find(strcmp(SACSCCmetrics{BOOTrep}.DCpeaks_legend,'DCpeak_CD'));
			SACSCCmetrics{BOOTrep}.CCCtfs_legend{3}='DCpeak_CD';
			% NOTE: peaks for ACFs are taken at 0 delay (by definition),
			% rather than at specified CCF delay
			SACSCCmetrics{BOOTrep}.CCCtfs(3) = SACSCCmetrics{BOOTrep}.DCpeak_AB(DC_CD_index) ...
				/ (sqrt(SACSCCmetrics{BOOTrep}.DCpeak_A(DC_0_index)*SACSCCmetrics{BOOTrep}.DCpeak_B(DC_0_index)));
			% 4+) User-passed list of delays to compute
			if isfield(paramsOUT,'UserDelays_usec')
				for i=1:length(paramsOUT.UserDelays_usec)
					DC_user_index=find(strcmp(SACSCCmetrics{BOOTrep}.DCpeaks_legend,sprintf('DCpeak_%d',round(paramsOUT.UserDelays_usec(i)))));
					SACSCCmetrics{BOOTrep}.CCCtfs_legend{end+1}=sprintf('DCpeak_%d',round(paramsOUT.UserDelays_usec(i)));
					% NOTE: peaks for ACFs are taken at 0 delay (by definition),
					% rather than at specified CCF delay
					SACSCCmetrics{BOOTrep}.CCCtfs(end+1) = SACSCCmetrics{BOOTrep}.DCpeak_AB(DC_user_index) ...
						/ (sqrt(SACSCCmetrics{BOOTrep}.DCpeak_A(DC_0_index)*SACSCCmetrics{BOOTrep}.DCpeak_B(DC_0_index)));
				end
			end

		else
			for i=length(SACSCCmetrics{BOOTrep}.DCpeaks_legend)
				SACSCCmetrics{BOOTrep}.CCCtfs(i) = NaN;
			end
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
		% using IFFT(0 delay) SCpeaks
		IFFTSC_0_index=find(strcmp(SACSCCmetrics{BOOTrep}.SCpeaks_legend,'IFFTraw_0'));
		SACSCCmetrics{BOOTrep}.CCCenvs_legend{end+1}='IFFTrawSC_0';
		SACSCCmetrics{BOOTrep}.CCCenvs(end+1) = (SACSCCmetrics{BOOTrep}.SCpeaks_AB(IFFTSC_0_index)-1)/ ...
			(sqrt((SACSCCmetrics{BOOTrep}.SCpeaks_A(IFFTSC_0_index)-1)*(SACSCCmetrics{BOOTrep}.SCpeaks_B(IFFTSC_0_index)-1)));
		% using IFFT(CD delay) SCpeaks
		IFFTSC_CD_index=find(strcmp(SACSCCmetrics{BOOTrep}.SCpeaks_legend,'IFFTraw_CD'));
		SACSCCmetrics{BOOTrep}.CCCenvs_legend{end+1}='IFFTrawSC_CD';
		% NOTE: peaks for ACFs are taken at 0 delay (by definition), rather
		% than at specified CCF delay 
		SACSCCmetrics{BOOTrep}.CCCenvs(end+1) = (SACSCCmetrics{BOOTrep}.SCpeaks_AB(IFFTSC_CD_index)-1)/ ...
			(sqrt((SACSCCmetrics{BOOTrep}.SCpeaks_A(IFFTSC_0_index)-1)*(SACSCCmetrics{BOOTrep}.SCpeaks_B(IFFTSC_0_index)-1)));
		% User-passed list of delays to compute
		if isfield(paramsOUT,'UserDelays_usec')
			for i=1:length(paramsOUT.UserDelays_usec)
				IFFTSC_user_index=find(strcmp(SACSCCmetrics{BOOTrep}.SCpeaks_legend,sprintf('IFFTraw_%d',round(paramsOUT.UserDelays_usec(i)))));
				SACSCCmetrics{BOOTrep}.CCCenvs_legend{end+1}=sprintf('IFFTrawSC_%d',round(paramsOUT.UserDelays_usec(i)));
				% NOTE: peaks for ACFs are taken at 0 delay (by definition),
				% rather than at specified CCF delay
				SACSCCmetrics{BOOTrep}.CCCenvs(end+1) = (SACSCCmetrics{BOOTrep}.SCpeaks_AB(IFFTSC_user_index)-1)/ ...
					(sqrt((SACSCCmetrics{BOOTrep}.SCpeaks_A(IFFTSC_0_index)-1)*(SACSCCmetrics{BOOTrep}.SCpeaks_B(IFFTSC_0_index)-1)));
			end
		end

	end

	SACSCCmetrics{BOOTrep}.NumDrivenSpikes=NumDrivenSpikes;
	SACSCCmetrics{BOOTrep}.AvgRate_sps=SACSCCmetrics_prePATCH100408{BOOTrep}.AvgRate_sps;


end

%% TALLY AVG
CCCtfsVECs=cell(size(SACSCCmetrics{1}.CCCtfs));
CCCtfsAVGs=zeros(size(SACSCCmetrics{1}.CCCtfs));
CCCtfsSTDs=zeros(size(SACSCCmetrics{1}.CCCtfs));
for i=1:length(CCCtfsVECs)
	CCCtfsVECs{i}=zeros(1,paramsOUT.BOOTSTRAP_Navgs);
	for BOOTrep4=1:paramsOUT.BOOTSTRAP_Navgs
		CCCtfsVECs{i}(BOOTrep4)=SACSCCmetrics{BOOTrep4}.CCCtfs(i);
	end
	CCCtfsAVGs(i)=mean(CCCtfsVECs{i});
	CCCtfsSTDs(i)=std(CCCtfsVECs{i});
end

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
SACSCCmetrics{paramsOUT.BOOTSTRAP_Navgs+1}.CCCtfsVECs=CCCtfsVECs;
SACSCCmetrics{paramsOUT.BOOTSTRAP_Navgs+1}.CCCtfsAVGs=CCCtfsAVGs;
SACSCCmetrics{paramsOUT.BOOTSTRAP_Navgs+1}.CCCtfsSTDs=CCCtfsSTDs;



plot_CCCanal_3(SACSCCfunctions,SACSCCmetrics,paramsOUT)
set(gcf,'Name',sprintf('Post-PATCH: %s',NEWfilename))

%% SAVE new CCCanal Data
disp(sprintf('SAVING updated (PATCHED) CCCanal_3 Data: %s.mat',NEWfilename))
eval(['save ' NEWfilename '.mat'])


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OUTPUT various metrics of interest for AVGed value
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SACpeakLIST={'SACpeak_max ','SACpeak_0   ','SACpeak_CD  '};
for i=1:length(SACpeakLIST)
	disp(sprintf('SACpeaks("%s") = A: %.2f;\tB: %.2f;\tAB: %.2f',SACpeakLIST{i}, ...
		SACSCCmetrics{end}.SACpeak_A(find(strcmp(deblank(SACpeakLIST{i}),SACSCCmetrics{end}.SACpeaks_legend))), ...
		SACSCCmetrics{end}.SACpeak_B(find(strcmp(deblank(SACpeakLIST{i}),SACSCCmetrics{end}.SACpeaks_legend))), ...
		SACSCCmetrics{end}.SCCpeak_AB(find(strcmp(deblank(SACpeakLIST{i}),SACSCCmetrics{end}.SACpeaks_legend))) ))
end
disp(sprintf(' '))
DCpeakLIST={'DCpeak_max ','DCpeak_0   ','DCpeak_CD  '};
for i=1:length(DCpeakLIST)
	disp(sprintf('DCpeaks("%s") = A: %.2f;\tB: %.2f;\tAB: %.2f',DCpeakLIST{i}, ...
		SACSCCmetrics{end}.DCpeak_A(find(strcmp(deblank(DCpeakLIST{i}),SACSCCmetrics{end}.DCpeaks_legend))), ...
		SACSCCmetrics{end}.DCpeak_B(find(strcmp(deblank(DCpeakLIST{i}),SACSCCmetrics{end}.DCpeaks_legend))), ...
		SACSCCmetrics{end}.DCpeak_AB(find(strcmp(deblank(DCpeakLIST{i}),SACSCCmetrics{end}.DCpeaks_legend))) ))
end
disp(sprintf(' '))
SCpeakLIST={'raw         ','adj: 0-300  ','adj: 0-CF   ','IFFTraw     ','IFFTraw_0   ','IFFTraw_CD  '};
for i=1:length(SCpeakLIST)
	disp(sprintf('SCpeaks("%s") = A: %.2f;\tB: %.2f;\tAB: %.2f',SCpeakLIST{i}, ...
		SACSCCmetrics{end}.SCpeaks_A(find(strcmp(deblank(SCpeakLIST{i}),SACSCCmetrics{end}.SCpeaks_legend))), ...
		SACSCCmetrics{end}.SCpeaks_B(find(strcmp(deblank(SCpeakLIST{i}),SACSCCmetrics{end}.SCpeaks_legend))), ...
		SACSCCmetrics{end}.SCpeaks_AB(find(strcmp(deblank(SCpeakLIST{i}),SACSCCmetrics{end}.SCpeaks_legend))) ))
end
disp(sprintf(' '))
CCCtfsLIST={'DCpeak_max ','DCpeak_0   ','DCpeak_CD  '};
for i=1:length(CCCtfsLIST)
	disp(sprintf('CCCtfs("%s") = %.2f',CCCtfsLIST{i}, ...
		SACSCCmetrics{end}.CCCtfs(find(strcmp(deblank(CCCtfsLIST{i}),SACSCCmetrics{end}.CCCtfs_legend)))))
end
disp(sprintf(' '))
CCCenvLIST={'rawSC','10-300, subBIAS','10-300, withBIAS','0-300, subBIAS','0-50, subBIAS','0-CF, subBIAS','IFFTrawSC','IFFTrawSC_0','IFFTrawSC_CD'};
for i=1:length(CCCenvLIST)
	disp(sprintf('CCCenv("%s") = %.2f',CCCenvLIST{i}, ...
		SACSCCmetrics{end}.CCCenvs(find(strcmp(CCCenvLIST{i},SACSCCmetrics{end}.CCCenvs_legend)))))
end

return;


