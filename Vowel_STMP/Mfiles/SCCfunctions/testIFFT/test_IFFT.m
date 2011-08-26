% FILE test_IFFT.m


%% Compute FFT 
FFT_A=fft((SACSCCfunctions.SUMCOR_A-1),paramsIN.Nfft_psd);


paramsIN.SACSCC_CF_Hz
[y,CF_index]=min(abs(SACSCCfunctions.freqVEC-paramsIN.SACSCC_CF_Hz));

FFT_A_adj=zeros(size(FFT_A));
FFT_A_adj(1:CF_index)=FFT_A(1:CF_index);
FFT_A_adj((length(FFT_A)-CF_index+1):end)=FFT_A((length(FFT_A)-CF_index+1):end);

adjSC=ifft(FFT_A_adj)+1;



figure(5)
subplot(211)
plot(SACSCCfunctions.freqVEC,abs(FFT_A))
hold on
plot(SACSCCfunctions.freqVEC,abs(FFT_A_adj),'r')

subplot(212)
% plot(SACSCCfunctions.delays_usec,SACSCCfunctions.SUMCOR_A) hold on
plot(SACSCCfunctions.SUMCOR_A) 
hold on
plot(real(adjSC),'r')


% plot(SACSCCfunctions.delays_usec,adjSC(1:length()




