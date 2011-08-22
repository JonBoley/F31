function [gain1,gain2,gain3,gain4,gain5] = OptimalGain2(directory,levels,impairment,strategy,phone,gains,note)
% [gain1,gain2,gain3,gain4,gain5] = OptimalGain2(directory,levels,impairment,strategy,phone,gains,note)
% gain1 optimizes short-term rate
% gain2 optimizes average rate
% gain3 optimizes rate (no windowing)
% gain4 optimizes envelope coding
% gain5 optimizes tfs coding
% example:
% OptimalGain(30:10:100,-40:10:40,'archive\65dBSPL\mixed\avg\phone12','Nov_24_08')

doplots=0;
w_spont = [.6 .2 .2]; %spont weights
existing_levels=[];
WarnMe = 0;

gain1 = zeros(length(levels),1);
gain2 = zeros(length(levels),1);
gain3 = zeros(length(levels),1);
gain4 = zeros(length(levels),1);
gain5 = zeros(length(levels),1);
gain1_range = zeros(length(levels),2);
gain2_range = zeros(length(levels),2);
gain3_range = zeros(length(levels),2);
gain4_range = zeros(length(levels),2);
gain5_range = zeros(length(levels),2);

level_index=1;
for level=levels
    if exist([directory num2str(level) 'dBSPL\' impairment '\' strategy '\phone' num2str(phone) '\0dBgain_' note '.mat'])
        existing_levels = [existing_levels, level];
        load([directory num2str(level) 'dBSPL\' impairment '\' strategy '\phone' num2str(phone) '\0dBgain_' note], '-regexp', '^neurogramB');
        index=1;
        for gain = gains
            if exist([directory num2str(level) 'dBSPL\' impairment '\' strategy '\phone' num2str(phone) '\' num2str(gain) 'dBgain_' note '.mat'])
                load([directory num2str(level) 'dBSPL\' impairment '\' strategy '\phone' num2str(phone) '\' num2str(gain) 'dBgain_' note], '-regexp', '^neurogramA','^Env','^Tfs');
                err1(level_index,index) = mean(mean(abs(neurogramB1-neurogramA1))); % short window
                err2(level_index,index) = mean(mean(abs(neurogramB2-neurogramA2))); % long window
                err3(level_index,index) = mean(mean(abs(neurogramB3-neurogramA3))); % no window
                
                err1(level_index,index) = err1(level_index,index)/mean(mean((neurogramB1+neurogramA1)/2)); %normalize by mean rate
                err2(level_index,index) = err2(level_index,index)/mean(mean((neurogramB2+neurogramA2)/2));
                err3(level_index,index) = err3(level_index,index)/mean(mean((neurogramB3+neurogramA3)/2));
                
                Env_temp = NaN*ones(size(Env,1),3); for i=1:size(Env,1), Env_temp(i,:)=Env{i,1}(:,1); end
                env(level_index,index) = mean(Env_temp*w_spont'); % take avg across CF (weighted SR's)
                BigEnv(index,:) = Env_temp*w_spont';
                %                 figure(98), %subplot(9,2,index), plot(Env);
                %                 hold on; plot(Env*w_spont'); hold off; %axis off;
                %                 title(['Env, ' num2str(gain) 'dB gain']);
                Tfs_temp = NaN*ones(size(Tfs,1),3); for i=1:size(Tfs,1), Tfs_temp(i,:)=Tfs{i,1}(:,1); end
                tfs(level_index,index) = mean(Tfs_temp*w_spont');
                BigTfs(index,:) = Tfs_temp*w_spont';
                %                 figure(99), %subplot(9,2,index), plot(Tfs);
                %                 title(['Tfs, ' num2str(gain) 'dB gain']);
                %                 hold on; plot(Tfs*w_spont'); hold off; %axis off;
                
                %         figure(1), subplot(length(gains),2,(2*index)-1), imagesc(neurogramB1'/max(max(neurogramB1)),[0 1]); axis off; title('Normal');
                %         subplot(length(gains),2,2*index), imagesc(neurogramA1'/max(max(neurogramB1)),[0 1]); axis off; title(sprintf('Gain = %ddB',gain));
                %         figure(2), subplot(length(gains),2,(2*index)-1), imagesc(neurogramB2'/max(max(neurogramB2)),[0 1]); axis off; title('Normal');
                %         subplot(length(gains),2,2*index), imagesc(neurogramA2'/max(max(neurogramB2)),[0 1]); axis off; title(sprintf('Gain = %ddB',gain));
                clear -regexp ^neurogramA ^Env ^Tfs
                index = index+1;
            else
                WarnMe = 1;  %display warning
            end % if gain exists
        end
        clear -regexp ^neurogramB
        [C,I] = min(err1(level_index,:));
        gain1(level_index) = gains(I);
        [r,c]=find(err1(level_index,:)<min(err1(level_index,:))+0.2*(max(err1(level_index,:))-min(err1(level_index,:)))); %top 20 percent
        gain1_range(level_index,:) = [min(gains(c)),max(gains(c))];
        [C,I] = min(err2(level_index,:));
        gain2(level_index) = gains(I);
        [r,c]=find(err2(level_index,:)<min(err2(level_index,:))+0.2*(max(err2(level_index,:))-min(err2(level_index,:))));
        gain2_range(level_index,:) = [min(gains(c)),max(gains(c))];
        [C,I] = min(err3(level_index,:));
        gain3(level_index) = gains(I);
        [r,c]=find(err3(level_index,:)<min(err3(level_index,:))+0.2*(max(err3(level_index,:))-min(err3(level_index,:))));
        gain3_range(level_index,:) = [min(gains(c)),max(gains(c))];
        [C,I] = max(env(level_index,:).*isfinite(env(level_index,:)));
        if ~isnan(C), gain4(level_index) = gains(I); end
        temp_rhoE = env(level_index,:).*isfinite(env(level_index,:));
        [r,c]=find(env(level_index,:)>min(temp_rhoE)+0.8*(max(temp_rhoE)-min(temp_rhoE)));
        if(c)
            gain4_range(level_index,:) = [min(gains(c)),max(gains(c))];
        else
            gain4_range(level_index,:) = [NaN,NaN];
        end
        [C,I] = max(tfs(level_index,:).*isfinite(tfs(level_index,:)));
        if ~isnan(C), gain5(level_index) = gains(I); end
        temp_rhoT = tfs(level_index,:).*isfinite(tfs(level_index,:));
        [r,c]=find(tfs(level_index,:)>min(temp_rhoT)+0.8*(max(temp_rhoT)-min(temp_rhoT)));
        if(c)
            gain5_range(level_index,:) = [min(gains(c)),max(gains(c))];
        else
            gain5_range(level_index,:) = [NaN,NaN];
        end
        level_index=level_index+1;
    end % if level exists
    
    gain1 = squeeze(gain1);
    gain2 = squeeze(gain2);
    gain3 = squeeze(gain3);
    gain4 = squeeze(gain4);
    gain5 = squeeze(gain5);
    gain1_range = squeeze(gain1_range);
    gain2_range = squeeze(gain2_range);
    gain3_range = squeeze(gain3_range);
    gain4_range = squeeze(gain4_range);
    gain5_range = squeeze(gain5_range);
    
    %     pause;
end

if(doplots)
    figure, subplot(2,2,1), errorbar(existing_levels,gain2,gain2-gain2_range(:,1),gain2_range(:,2)-gain2);
    axis([0 120 min(gains(1),gains(end)) max(gains(1),gains(end))]);
    title(sprintf('NAL Gain Adjustment\n(Average Discharge Rate)'));
    xlabel('Input SPL (dB)'); ylabel('Gain Adjustment (dB)');
    
    subplot(2,2,2), errorbar(existing_levels,gain1,gain1-gain1_range(:,1),gain1_range(:,2)-gain1);
    axis([0 120 min(gains(1),gains(end)) max(gains(1),gains(end))]);
    title(sprintf('NAL Gain Adjustment\n(Short-Term Discharge Rate)'));
    xlabel('Input SPL (dB)'); ylabel('Gain Adjustment (dB)');
    
    % subplot(2,3,3), plot(existing_levels,gain3,'x','MarkerSize',8);
    % axis([0 120 min(gains(1),gains(end)) max(gains(1),gains(end))]);
    % title(sprintf('NAL Gain Adjustment\n(Expected Discharge Rate)'));
    % xlabel('Input SPL (dB)'); ylabel('Gain Adjustment (dB)');
    
    subplot(2,2,3), errorbar(existing_levels,gain4,gain4-gain4_range(:,1),gain4_range(:,2)-gain4);
    axis([0 120 min(gains(1),gains(end)) max(gains(1),gains(end))]);
    title(sprintf('NAL Gain Adjustment\n(Envelope Coding)'));
    xlabel('Input SPL (dB)'); ylabel('Gain Adjustment (dB)');
    
    subplot(2,2,4), errorbar(existing_levels,gain5,gain5-gain5_range(:,1),gain5_range(:,2)-gain5);
    axis([0 120 min(gains(1),gains(end)) max(gains(1),gains(end))]);
    title(sprintf('NAL Gain Adjustment\n(TFS Coding)'));
    xlabel('Input SPL (dB)'); ylabel('Gain Adjustment (dB)');
end

if WarnMe, warndlg('Some gain values did not exist!'); end

