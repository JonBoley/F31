%% This version plots normal, impaired, and aided all on the same plot

% dbstop if warning MATLAB:CDwarn

DoSNR = 1;
if DoSNR
    F1strings = {...
        '12-May-2010_snr_f1_nh.mat',...
        '12-May-2010_snr_f1_snhl.mat',...
        '19-May-2010_snr_f1_snhl_aided.mat',...
        '01-Jun-2010_snr_f1_snhl_dsl.mat'};
    F2strings = {...
        '19-May-2010_snr_f2_nh.mat',...
        '19-May-2010_snr_f2_snhl.mat',...
        '25-May-2010_snr_f2_snhl_aided.mat',...
        '01-Jun-2010_snr_f2_snhl_dsl.mat'};
    F3strings = {...
        '24-May-2010_snr_f3_nh.mat',...
        '24-May-2010_snr_f3_snhl.mat',...
        '24-May-2010_snr_f3_snhl_aided.mat',...
        '01-Jun-2010_snr_f3_snhl_dsl.mat'};
else
    F1strings = {...
        '28-May-2010_levels_f1_nh.mat',...
        '28-May-2010_levels_f1_snhl.mat',...
        '28-May-2010_levels_f1_snhl_aided.mat',...
        '01-Jun-2010_levels_f1_snhl_dsl.mat'};
    F2strings = {...
        '28-May-2010_levels_f2_nh.mat',...
        '28-May-2010_levels_f2_snhl.mat',...
        '28-May-2010_levels_f2_snhl_aided.mat',...
        '01-Jun-2010_levels_f2_snhl_dsl.mat'};
    F3strings = {...
        '28-May-2010_levels_f3_nh.mat',...
        '28-May-2010_levels_f3_snhl.mat',...
        '28-May-2010_levels_f3_snhl_aided.mat',...
        '01-Jun-2010_levels_f3_snhl_dsl.mat'};
end

for hearing = 1:4 % {'nh','snhl','snhl_aided','nh_aided'}

    hearingstrings = {'Normal','SNHL','SNHL + Linear Gain','SNHL + Compression'};
    %% F1
    load(F1strings{hearing});
    lines = {'go-','rs-','b^-','bv:'};
    LocalRho = squeeze(Rho_vDeltaCF);
    LocalCD = squeeze(CD_vDeltaCF);
    LocalCD = LocalCD*1e-6/(1/(CF_kHz(1)*1e3)); %convert to CF cycles
    for i=1:size(LocalRho,1)
        figure(1000-i), subplot(2,1,1), hold on;
        set(gcf,'WindowStyle','docked');
        LocalRho(i,midCF_index)=NaN;
        h_legend(hearing,i)=semilogx([CF_kHz(2:midCF_index) CF_kHz(1) CF_kHz(midCF_index+1:end)]',LocalRho(i,:),lines{hearing},'LineWidth',3);
        legend(h_legend(:,i),hearingstrings{1:hearing});
        if hearing==1
            plot([CF_kHz(1) CF_kHz(1)],[min(min(Rho_vDeltaCF)) min(1,max(max(Rho_vDeltaCF)))],'k:');
            text(CF_kHz(1),0.05+min(1,max(max(Rho_vDeltaCF))),FeaturesText{featureNum},...
                'HorizontalAlignment','Center');
        end
    end
    tempCD_vDeltaCF=cell(size(LocalCD,1),1);
    tempCF_kHz=cell(size(LocalCD,1),1);
    freqs = [CF_kHz(2:midCF_index) CF_kHz(1) CF_kHz(midCF_index+1:end)]';
    for i=1:size(LocalCD,1)
        figure(1000-i), subplot(2,1,2), hold on;
        h = semilogx(freqs,LocalCD(i,:),lines{hearing},'LineWidth',3); hold on;
        reply = input('Does the CD look okay? Y/N [Y]: ', 's');
        if isempty(reply), reply = 'Y'; end
        if ~strcmp(reply,'Y')
            reply = input('Which indices should be ignored?');
            tempCD_vDeltaCF{i}=LocalCD(i,setdiff(1:length(LocalCD(i,:)),reply))';
            tempCF_kHz{i}=freqs(setdiff(1:length(LocalCD(i,:)),reply));
            delete(h);
            semilogx(tempCF_kHz{i},squeeze(tempCD_vDeltaCF{i}),lines{hearing},'LineWidth',3);
            tempCD_vDeltaCF{i}=LocalCD(i,...
                setdiff((fix(length(LocalCD(i,:))/2)):(fix(length(LocalCD(i,:))/2)+2),reply))';
            tempCF_kHz{i}=freqs(...
                setdiff((fix(length(LocalCD(i,:))/2)):(fix(length(LocalCD(i,:))/2)+2),reply));
%             warning('Matlab:CDwarn', 'Please edit CD');
        else
            % just the middle 3 points please
            tempCD_vDeltaCF{i}=LocalCD(i,...
                (fix(length(LocalCD(i,:))/2)):(fix(length(LocalCD(i,:))/2)+2))';
            tempCF_kHz{i}=freqs((fix(length(freqs)/2)):(fix(length(freqs)/2)+2)); % just the middle 3 points
        end
    end
    text(CF_kHz(1),500,FeaturesText{featureNum},...
        'HorizontalAlignment','Center');
    
    figure(1000), subplot(2,3,1), hold on;
    set(gcf,'WindowStyle','docked');
    for i=1:size(LocalRho,1) % for each SNR
        rho_ref = mean([LocalRho(i,midCF_index-1),...
            LocalRho(i,midCF_index+1)]);
        width_index2 = interp1(LocalRho(i,midCF_index+1:end),CF_kHz(midCF_index+1:end),0.80*rho_ref);
        width_index1 = interp1(LocalRho(i,1:midCF_index-1),CF_kHz(2:midCF_index),0.80*rho_ref);
        width_oct(i) = log2(width_index2/width_index1);
        if isnan(width_oct(i)) %clip measurement
            width_oct(i)=1.01;
        end
    end
    if DoSNR
        Xnumbers = [12 SNRs(2:end)];
        plot([12 SNRs(2:end)],width_oct,lines{hearing},'LineWidth',3);
        for i=1:length(width_oct)
           if width_oct(i)==1.01
               text(Xnumbers(i),1,'\uparrow','FontWeight','bold','VerticalAlignment','bottom');
           end
        end
        set(gca,'XTick',[fliplr(SNRs(2:end)) 12]);
        set(gca,'XTickLabel','-6|0|6|Quiet');
        set(gca,'Xdir','reverse');
    else
        plot(Levels,width_oct,lines{hearing},'LineWidth',3);
    end
    ylabel('\rho width (octaves)');
    axis([-6 18 0.2 1.01]);
    title('F1');
    
    figure(1000), subplot(2,3,4), hold on;
    for i=1:size(LocalRho,1) % for each SNR
        if isempty(tempCD_vDeltaCF{i})
            Cfit=polyfit(log2([CF_kHz(2:midCF_index) CF_kHz(1) CF_kHz(midCF_index+1:end)]'/CF_kHz(1)),...
                LocalCD(i,:),1); % Cfit=[slope, intercept]
        else
            Cfit=polyfit(log2(tempCF_kHz{i}/CF_kHz(1)),...
                tempCD_vDeltaCF{i},1); % Cfit=[slope, intercept]
            if length(tempCD_vDeltaCF{i})<3 % if we don't have all 3 center points
                Cfit(1)=NaN;
            end
        end
        slope(i) = Cfit(1);
    end
    if DoSNR
        plot([12 SNRs(2:end)],-slope,lines{hearing},'LineWidth',3);
        set(gca,'XTick',[fliplr(SNRs(2:end)) 12]);
        set(gca,'XTickLabel','-6|0|6|Quiet');
        set(gca,'Xdir','reverse');
        xlabel('SNR (dB)');
    else 
        plot(Levels,-slope,lines{hearing},'LineWidth',3);
        xlabel('Level (dB SPL)');
    end
    ylabel('-CD slope (cycles/octave)');
	axis([-6 18 0 5]);
    
    clearvars -except hearing h_legend F1strings F2strings F3strings DoSNR;
    
    %% F2
    load(F2strings{hearing});
    lines = {'go-','rs-','b^-','bv:'};
    LocalRho = squeeze(Rho_vDeltaCF);
    LocalCD = squeeze(CD_vDeltaCF);
    LocalCD = LocalCD*1e-6/(1/(CF_kHz(1)*1e3)); %convert to CF cycles
    for i=1:size(LocalRho,1)
        figure(1000-i), subplot(2,1,1), hold on;
        set(gcf,'WindowStyle','docked');
        LocalRho(i,midCF_index)=NaN;
        semilogx([CF_kHz(2:midCF_index) CF_kHz(1) CF_kHz(midCF_index+1:end)]',LocalRho(i,:),lines{hearing},'LineWidth',3);
        if hearing==1
            plot([CF_kHz(1) CF_kHz(1)],[min(min(Rho_vDeltaCF)) min(1,max(max(Rho_vDeltaCF)))],'k:');
            text(CF_kHz(1),0.05+min(1,max(max(Rho_vDeltaCF))),FeaturesText{featureNum},...
                'HorizontalAlignment','Center');
        end
    end
    tempCD_vDeltaCF=cell(size(LocalCD,1),1);
    tempCF_kHz=cell(size(LocalCD,1),1);
    freqs = [CF_kHz(2:midCF_index) CF_kHz(1) CF_kHz(midCF_index+1:end)]';
    for i=1:size(LocalCD,1)
        figure(1000-i), subplot(2,1,2), hold on;
        h = semilogx(freqs,LocalCD(i,:),lines{hearing},'LineWidth',3); hold on;
        reply = input('Does the CD look okay? Y/N [Y]: ', 's');
        if isempty(reply), reply = 'Y'; end
        if ~strcmp(reply,'Y')
            reply = input('Which indices should be ignored?');
            tempCD_vDeltaCF{i}=LocalCD(i,setdiff(1:length(LocalCD(i,:)),reply))';
            tempCF_kHz{i}=freqs(setdiff(1:length(LocalCD(i,:)),reply));
            delete(h);
            semilogx(tempCF_kHz{i},squeeze(tempCD_vDeltaCF{i}),lines{hearing},'LineWidth',3);
            tempCD_vDeltaCF{i}=LocalCD(i,...
                setdiff((fix(length(LocalCD(i,:))/2)):(fix(length(LocalCD(i,:))/2)+2),reply))';
            tempCF_kHz{i}=freqs(...
                setdiff((fix(length(LocalCD(i,:))/2)):(fix(length(LocalCD(i,:))/2)+2),reply));
%             warning('Matlab:CDwarn', 'Please edit CD');
        else
            % just the middle 3 points please
            tempCD_vDeltaCF{i}=LocalCD(i,...
                (fix(length(LocalCD(i,:))/2)):(fix(length(LocalCD(i,:))/2)+2))';
            tempCF_kHz{i}=freqs((fix(length(freqs)/2)):(fix(length(freqs)/2)+2)); % just the middle 3 points
        end
    end
    text(CF_kHz(1),500,FeaturesText{featureNum},...
        'HorizontalAlignment','Center');
    
    figure(1000), subplot(2,3,2), hold on;
    set(gcf,'WindowStyle','docked');
    for i=1:size(LocalRho,1) % for each SNR
        rho_ref = mean([LocalRho(i,midCF_index-1),...
            LocalRho(i,midCF_index+1)]);
        width_index2 = interp1(LocalRho(i,midCF_index+1:end),CF_kHz(midCF_index+1:end),0.80*rho_ref);
        width_index1 = interp1(LocalRho(i,1:midCF_index-1),CF_kHz(2:midCF_index),0.80*rho_ref);
        width_oct(i) = log2(width_index2/width_index1);
        if isnan(width_oct(i)) %clip measurement
            width_oct(i)=1.01;
        end
    end
    if DoSNR
        Xnumbers = [12 SNRs(2:end)];
        plot([12 SNRs(2:end)],width_oct,lines{hearing},'LineWidth',3);
        for i=1:length(width_oct)
           if width_oct(i)==1.01
               text(Xnumbers(i),1,'\uparrow','FontWeight','bold','VerticalAlignment','bottom');
           end
        end
        set(gca,'XTick',[fliplr(SNRs(2:end)) 12]);
        set(gca,'XTickLabel','-6|0|6|Quiet');
        set(gca,'Xdir','reverse');
    else
        plot(Levels,width_oct,lines{hearing},'LineWidth',3);
    end
%     ylabel('\rho width (octaves)');
    axis([-6 18 0.2 1.01]);
    title('F2');
    
    figure(1000), subplot(2,3,5), hold on;
    for i=1:size(LocalRho,1) % for each SNR/Level
        if isempty(tempCD_vDeltaCF{i})
            Cfit=polyfit(log2([CF_kHz(2:midCF_index) CF_kHz(1) CF_kHz(midCF_index+1:end)]'/CF_kHz(1)),...
                LocalCD(i,:),1); % Cfit=[slope, intercept]
        else
            Cfit=polyfit(log2(tempCF_kHz{i}/CF_kHz(1)),...
                tempCD_vDeltaCF{i},1); % Cfit=[slope, intercept]
            if length(tempCD_vDeltaCF{i})<3 % if we don't have all 3 center points
                Cfit(1)=NaN;
            end
        end
        slope(i) = Cfit(1);
    end
    if DoSNR
        plot([12 SNRs(2:end)],-slope,lines{hearing},'LineWidth',3);
        set(gca,'XTick',[fliplr(SNRs(2:end)) 12]);
        set(gca,'XTickLabel','-6|0|6|Quiet');
        set(gca,'Xdir','reverse');
        xlabel('SNR (dB)');
    else 
        plot(Levels,-slope,lines{hearing},'LineWidth',3);
        xlabel('Level (dB SPL)');
    end
%     ylabel('-CD slope (cycles/octave)');
	axis([-6 18 0 5]);
    
    clearvars -except hearing h_legend F1strings F2strings F3strings DoSNR;
    
    %% F3
    load(F3strings{hearing});
    lines = {'go-','rs-','b^-','bv:'};
    LocalRho = squeeze(Rho_vDeltaCF);
    LocalCD = squeeze(CD_vDeltaCF);
    LocalCD = LocalCD*1e-6/(1/(CF_kHz(1)*1e3)); %convert to CF cycles
    for i=1:size(LocalRho,1)
        figure(1000-i), subplot(2,1,1), hold on;
        set(gcf,'WindowStyle','docked');
        LocalRho(i,midCF_index)=NaN;
        semilogx([CF_kHz(2:midCF_index) CF_kHz(1) CF_kHz(midCF_index+1:end)]',LocalRho(i,:),lines{hearing},'LineWidth',3);
        if hearing==1
            plot([CF_kHz(1) CF_kHz(1)],[min(min(Rho_vDeltaCF)) min(1,max(max(Rho_vDeltaCF)))],'k:');
            text(CF_kHz(1),0.05+min(1,max(max(Rho_vDeltaCF))),FeaturesText{featureNum},...
                'HorizontalAlignment','Center');
        end
    end
    tempCD_vDeltaCF=cell(size(LocalCD,1),1);
    tempCF_kHz=cell(size(LocalCD,1),1);
    freqs = [CF_kHz(2:midCF_index) CF_kHz(1) CF_kHz(midCF_index+1:end)]';
    for i=1:size(LocalCD,1)
        figure(1000-i), subplot(2,1,2), hold on;
        h = semilogx(freqs,LocalCD(i,:),lines{hearing},'LineWidth',3); hold on;
        reply = input('Does the CD look okay? Y/N [Y]: ', 's');
        if isempty(reply), reply = 'Y'; end
        if ~strcmp(reply,'Y')
            reply = input('Which indices should be ignored?');
            tempCD_vDeltaCF{i}=LocalCD(i,setdiff(1:length(LocalCD(i,:)),reply))';
            tempCF_kHz{i}=freqs(setdiff(1:length(LocalCD(i,:)),reply));
            delete(h);
            semilogx(tempCF_kHz{i},squeeze(tempCD_vDeltaCF{i}),lines{hearing},'LineWidth',3);
            tempCD_vDeltaCF{i}=LocalCD(i,...
                setdiff((fix(length(LocalCD(i,:))/2)):(fix(length(LocalCD(i,:))/2)+2),reply))';
            tempCF_kHz{i}=freqs(...
                setdiff((fix(length(LocalCD(i,:))/2)):(fix(length(LocalCD(i,:))/2)+2),reply));
%             warning('Matlab:CDwarn', 'Please edit CD');
        else
            % just the middle 3 points please
            tempCD_vDeltaCF{i}=LocalCD(i,...
                (fix(length(LocalCD(i,:))/2)):(fix(length(LocalCD(i,:))/2)+2))';
            tempCF_kHz{i}=freqs((fix(length(freqs)/2)):(fix(length(freqs)/2)+2)); % just the middle 3 points
        end
    end
    text(CF_kHz(1),500,FeaturesText{featureNum},...
        'HorizontalAlignment','Center');
    
    for i=1:size(LocalRho,1)
        figure(1000-i), subplot(2,1,1), hold off;
        if DoSNR
            title(num2str(SNRs(i)));
        else
            title(num2str(Levels(i)));
        end
        
%         legend('Normal','SNHL','SNHL + Linear Gain','SNHL + Compression');
        xlabel('CF (kHz)');
        ylabel(['\rho (re F_x)']);
        axis([0 4 0.4 1.2]);
        set(gca,'XScale','log');
        subplot(2,1,2), hold off;
        xlabel('CF (kHz)');
        ylabel(['Characteristic Delay (cycles re F_x)']);
        axis([0 4 -1.5 1.5]);
        set(gca,'XScale','log');
    end
    
    figure(1000), subplot(2,3,3), hold on;
    for i=1:size(LocalRho,1) % for each SNR
        rho_ref = mean([LocalRho(i,midCF_index-1),...
            LocalRho(i,midCF_index+1)]);
        width_index2 = interp1(LocalRho(i,midCF_index+1:end),CF_kHz(midCF_index+1:end),0.80*rho_ref);
        width_index1 = interp1(LocalRho(i,1:midCF_index-1),CF_kHz(2:midCF_index),0.80*rho_ref);
        width_oct(i) = log2(width_index2/width_index1);
        if isnan(width_oct(i)) %clip measurement
            width_oct(i)=1.01;
        end
    end
    if DoSNR
        Xnumbers = [12 SNRs(2:end)];
        plot([12 SNRs(2:end)],width_oct,lines{hearing},'LineWidth',3);
        for i=1:length(width_oct)
           if width_oct(i)==1.01
               text(Xnumbers(i),1,'\uparrow','FontWeight','bold','VerticalAlignment','bottom');
           end
        end
        set(gca,'XTick',[fliplr(SNRs(2:end)) 12]);
        set(gca,'XTickLabel','-6|0|6|Quiet');
        set(gca,'Xdir','reverse');
    else
        plot(Levels,width_oct,lines{hearing},'LineWidth',3);
    end
    axis([-6 18 0.2 1.01]);
    legend('NH','SNHL','Linear','NonLin');
    title('F3');
    
    figure(1000), subplot(2,3,6), hold on;
    for i=1:size(LocalCD,1) % for each SNR
        if isempty(tempCD_vDeltaCF{i})
            Cfit=polyfit(log2([CF_kHz(2:midCF_index) CF_kHz(1) CF_kHz(midCF_index+1:end)]'/CF_kHz(1)),...
                LocalCD(i,:),1); % Cfit=[slope, intercept]
        else
            Cfit=polyfit(log2(tempCF_kHz{i}/CF_kHz(1)),...
                tempCD_vDeltaCF{i},1); % Cfit=[slope, intercept]
            if length(tempCD_vDeltaCF{i})<3 % if we don't have all 3 center points
                Cfit(1)=NaN;
            end
        end
        slope(i) = Cfit(1);
    end
    if DoSNR
        plot([12 SNRs(2:end)],-slope,lines{hearing},'LineWidth',3);
        set(gca,'XTick',[fliplr(SNRs(2:end)) 12]);
        set(gca,'XTickLabel','-6|0|6|Quiet');
        set(gca,'Xdir','reverse');
        xlabel('SNR (dB)');
    else
        plot(Levels,-slope,lines{hearing},'LineWidth',3);
        xlabel('Level (dB SPL)');
    end
	axis([-6 18 0 5]);
    
    clearvars -except hearing h_legend F1strings F2strings F3strings DoSNR;
    
end