%%
% assume this figure exludes any non-normal data
% (e.g., ignore any experiment for which >25% of Q10 values are
% outside of the 5-95% range in Kale & Heinz (2010) )
figFileName = 'C:\Users\Jon\Google Drive\Dissertation\Aim1\ANthresholds_inliers.fig';

open(figFileName);
hFig = gcf;

axesObjs = get(hFig, 'Children'); %axes handles
dataObjs = get(axesObjs(2), 'Children');

normalThresh = [];
impairedThresh = [];
for ii=1:numel(dataObjs)
    if strcmp(get(dataObjs(ii),'Marker'),'.') && ...
            all(abs(get(dataObjs(ii),'Color')-[0.5 0.5 0.5])<0.01)
        % gray (normal data)
        normalThresh = [normalThresh ...
            [get(dataObjs(ii),'XData');get(dataObjs(ii),'YData')]];
    elseif strcmp(get(dataObjs(ii),'Marker'),'.') && ...
            all(abs(get(dataObjs(ii),'Color')-[1 0 0])<0.01)
        % red (impaired)
        impairedThresh = [impairedThresh ...
            [get(dataObjs(ii),'XData');get(dataObjs(ii),'YData')]];
    end
end

figure, subplot(211),
plot(normalThresh(1,:),normalThresh(2,:),'k.');
hold on; plot(impairedThresh(1,:),impairedThresh(2,:),'r.');
set(gca,'XScale','log'); 
set(gca,'XTick',[0.25 0.5 1 2 4 8 10]);
set(gca,'FontName','Arial','FontSize',14);
xmin=0.25; xmax=10; ymin=-20; ymax=110;
axis([xmin xmax ymin ymax]);
ylabel('Threshold (dB SPL)','FontName','Arial','FontSize',14); 
xlabel('Frequency (kHz)','FontName','Arial','FontSize',14);

%% plot trends
freqs_kHz_std=[250,500,1000,2000,3000,4000,6000]*1e-3;
MinThresh = NaN*ones(length(freqs_kHz_std),1);
AvgThresh = NaN*ones(length(freqs_kHz_std),1);
MinThresh_Normal = NaN*ones(length(freqs_kHz_std),1);
AvgThresh_Normal = NaN*ones(length(freqs_kHz_std),1);
bw = 1/4; %octaves

for ii=1:numel(freqs_kHz_std)
    localThreshes = normalThresh(2,...
        (normalThresh(1,:)>=freqs_kHz_std(ii)*2^(-bw) &...
        normalThresh(1,:)<=freqs_kHz_std(ii)*2^(bw)));
    if ~isempty(localThreshes), AvgThresh_Normal(ii)=mean(localThreshes); end
    if ~isempty(localThreshes), MinThresh_Normal(ii)=min(localThreshes); end
    
    localThreshes = impairedThresh(2,...
        (impairedThresh(1,:)>=freqs_kHz_std(ii)*2^(-bw) &...
        impairedThresh(1,:)<=freqs_kHz_std(ii)*2^(bw)));
    if ~isempty(localThreshes), AvgThresh(ii)=mean(localThreshes); end
    if ~isempty(localThreshes), MinThresh(ii)=min(localThreshes); end
end

plot(freqs_kHz_std,AvgThresh_Normal,'ko:','Linewidth',2);
plot(freqs_kHz_std,MinThresh_Normal,'ko-','Linewidth',2);
plot(freqs_kHz_std,AvgThresh,'ro:','Linewidth',2);
plot(freqs_kHz_std,MinThresh,'ro-','Linewidth',2);


%% plot Q10
dataObjs = get(axesObjs(1), 'Children');

normalQ = [];
impairedQ = [];
for ii=1:numel(dataObjs)
    if strcmp(get(dataObjs(ii),'Marker'),'.') && ...
            all(abs(get(dataObjs(ii),'Color')-[0.5 0.5 0.5])<0.01)
        % gray (normal data)
        normalQ = [normalQ ...
            [get(dataObjs(ii),'XData');get(dataObjs(ii),'YData')]];
    elseif strcmp(get(dataObjs(ii),'Marker'),'.') && ...
            all(abs(get(dataObjs(ii),'Color')-[1 0 0])<0.01)
        % red (impaired)
        impairedQ = [impairedQ ...
            [get(dataObjs(ii),'XData');get(dataObjs(ii),'YData')]];
    end
end

subplot(212), plot(normalQ(1,:),normalQ(2,:),'k.');
hold on; plot(impairedQ(1,:),impairedQ(2,:),'r.');
set(gca,'XScale','log','YScale','log'); 
set(gca,'XTick',[0.25 0.5 1 2 4 8 10]);
set(gca,'YTick',[1 10],'YTickLabel',[1 10]);
set(gca,'FontName','Arial','FontSize',14);
xmin=0.25; xmax=10; ymin=0.5; ymax=10;
axis([xmin xmax ymin ymax]);
ylabel('Q10','FontName','Arial','FontSize',14); 
xlabel('Frequency (kHz)','FontName','Arial','FontSize',14);

load normema;
%%%%%%% Add lines from Kale and Heinz (2010) normal-hearing CHIN data
QlowM97(:,1)=[.4 10]';  QhighM97(:,1)=QlowM97(:,1);
% 5th and 95th percentiles from Kale and Heinz (2009)
QlowM97(:,2)=10^.1572*QlowM97(:,1).^.1416;
QhighM97(:,2)=10^.7989*QlowM97(:,1).^.1416;
plot(QlowM97(:,1),QlowM97(:,2),'k-','Linewidth',2);
plot(QhighM97(:,1),QhighM97(:,2),'k-','Linewidth',2);


%% Calculate audiogram
f = freqs_kHz_std; %kHz
audiogram = AvgThresh_Normal - AvgThresh;
audiogram2 = MinThresh_Normal - MinThresh;
figure, semilogx(freqs_kHz_std,audiogram,'bo:'); hold on;
semilogx(freqs_kHz_std,audiogram2,'bo-');
legend('Avg','Min');
xlim([min(freqs_kHz_std) max(freqs_kHz_std)]); ylim([-100 0]);
set(gca,'XTick',freqs_kHz_std,'XTickLabel',freqs_kHz_std);
title('Audiogram'); ylabel('dB HL'); xlabel('Frequency (kHz)');

% ABR
semilogx([500 1000 2000 4000 8000]/1e3,[-15.5740 -23.7651 -19.7825 -16.2083 -9.4868],'go-');

% audiogram
audiogram = [16   18   20    9    9]; % based on avg animal data (500OBN exposure)
freqs_Hz = [500,1000,2000,4000,6000];
semilogx(freqs_Hz/1e3,-audiogram,'k-');

