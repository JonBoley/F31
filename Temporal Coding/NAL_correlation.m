%% OptimalGain(70,20:-5:-20,'archive\phone12\','mixed','Dec_02_08');

% flip things around so that a tie goes to the minimal gain
[c,i]=max(fliplr(min(1,BigEnv).*isfinite(BigEnv)));
OpGains = gains(length(gains)-i+1);

FREQUENCIES = [250,500,1000,2000,3000,4000,6000]; % Hz
dBLoss = [20 20 30 40 45 50 50]; % Mild Hearing Loss (according to Bruce (ISAAR 2007)
k_NAL = [-18 -8 1 -1 -2 -2 -2]; % dB
H_3FA = sum(dBLoss(2:4)) / 3; % sum up loss at 500Hz,1kHz,2kHz
X = 0.15*H_3FA;
R = 0.31; % NAL-RP formula  (not quite half-gain rule)
NAL_IG = X + R.*dBLoss + k_NAL; % insertion gain
NAL_IG = max(0,NAL_IG); % no negative gain
lowCF=250; %Hz
highCF=4000; %Hz
numCFs=30;
spread=log2(highCF/lowCF); %number of octaves to span
CF_kHz=(250*2.^(0:spread/(numCFs-1):spread))/1000;
NAL = interp1(FREQUENCIES,NAL_IG,CF_kHz(1:30)*1000);

figure,scatter(OpGains,NAL);
xlabel('optimal gain adjustment'); ylabel('prescriptive gain');

[RHO,PVAL] = corrcoef(OpGains,NAL);
legend('Mixed Damage')
title('Optimal ENV coding at each frequency (70dBSPL overall)');

disp(sprintf('\rho=%1.2f',RHO(1,2)));
disp(sprintf('p=%1.5f',PVAL(1,2)));
