%%
W_window = .000256*100000;
w_spont = [.6 .2 .2]; %[60%HSR 20%MSR 20%LSR]
for level_index = 1:8
    for gain_index = 1:17
        neurogramA = zeros(length(SynOutA{1,level_index,gain_index,1,1}),30);%impaired
        neurogramB = zeros(length(SynOutB{1,level_index,gain_index,1,1}),30);%normal
        for i=1:30
            for j=1:3 % for SR's
                if ~isempty(SynOutA{1,level_index,gain_index,i,j})
                    neurogramA(:,i) = neurogramA(:,i) + w_spont(j)*SynOutA{1,level_index,gain_index,i,j}';
                    neurogramB(:,i) = neurogramB(:,i) + w_spont(j)*SynOutB{1,level_index,gain_index,i,j}';
                end
            end
            neurogramA(:,i) = filter(hann(W_window),1,neurogramA(:,i));
            neurogramB(:,i) = filter(hann(W_window),1,neurogramB(:,i));
        end
        Mydiff(level_index,gain_index) = mean(mean(abs(neurogramB-neurogramA)));
    end
end
plot((min(Mydiff,[],1)))
% OptGain = gains(min(Mydiff,[],1));
% figure;imagesc((neurogram)');
