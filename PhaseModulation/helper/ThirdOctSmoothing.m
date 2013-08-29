function spectrum_smoothed = ThirdOctSmoothing(spectrum_orig,freqs)
% Applies 1/3-Octave smoothing to spectral coefficients

N=1/3; % smoothing factor (octaves)

isColMatrix = find(size(spectrum_orig)==min(size(spectrum_orig)))-1;
if ~isColMatrix
    spectrum_orig=spectrum_orig';
    freqs=freqs';
end
    
for m=1:size(spectrum_orig,2)
    for i=1:size(spectrum_orig,1)
        minIndex = max(1,find(freqs<=freqs(i)*2^-N,1,'last'));
        if isempty(minIndex), minIndex=1; end
        maxIndex = min(length(spectrum_orig(:,m)),find(freqs>=freqs(i)*2^N,1,'first'));
        if isempty(maxIndex), maxIndex=length(spectrum_orig(:,m)); end
        
        spectrum_smoothed(i,m) = mean(spectrum_orig(minIndex:maxIndex,m));
    end
end

if ~isColMatrix
    spectrum_smoothed=spectrum_smoothed';
end
