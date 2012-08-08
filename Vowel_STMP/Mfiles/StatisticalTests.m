%%
CFrange = [0.1 10];%4*[2^(-0.5) 1*2^(0.5)]; %kHz
FeatureNum = [1 2];
SNRindex = [1 2 3];

NHindex = find(nhCFs>=CFrange(1) & nhCFs<=CFrange(end));
SNHLindex = find(snhlCFs>=CFrange(1) & snhlCFs<=CFrange(end));

% matrix{1}=[]; matrix{2}=[];
% for m=FeatureNum
%     for n=SNRindex
%         matrix{1} = [NHmatrix array_cd.nh{m,n}(NHindex)];
%         matrix{2} = [SNHLmatrix array_cd.snhl{m,n}(SNHLindex)];
%     end
% end

matrix{1}=[]; matrix{2}=[]; matrix{3}=[];
for m=FeatureNum
    for n=SNRindex
        matrix{n} = [matrix{n} array_cd.snhl{m,n}(SNHLindex)];
    end
end

[p_rs,h_rs,stats_rs] = ranksum(matrix{1}(~isnan(matrix{1})),matrix{2}(~isnan(matrix{2})));

if h_rs
    fprintf('Significant difference found! (p=%1.5f)\n',p_rs);
else
    fprintf('No significant difference found (p=%1.5f)\n',p_rs);
end
