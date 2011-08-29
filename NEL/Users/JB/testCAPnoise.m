%File TESTNOISE
%

x2=a0002_ABR_1998;
x3=a0003_ABR_1998;
x4=a0004_ABR_1998;
x5=a0005_ABR_1998;
x6=a0006_ABR_1998;
x7=a0007_ABR_1998;


RMS=sqrt(mean(x2.AD_Data.AD_Avg_V.^2));
% figure(2)
% plot(x2.AD_Data.AD_Avg_V)
% hold on
NoiseFACT=1.25;
x2.AD_Data.AD_Avg_V=x2.AD_Data.AD_Avg_V+NoiseFACT*RMS*randn(size(x2.AD_Data.AD_Avg_V));
x3.AD_Data.AD_Avg_V=x3.AD_Data.AD_Avg_V+NoiseFACT*RMS*randn(size(x2.AD_Data.AD_Avg_V));
x4.AD_Data.AD_Avg_V=x4.AD_Data.AD_Avg_V+NoiseFACT*RMS*randn(size(x2.AD_Data.AD_Avg_V));
x5.AD_Data.AD_Avg_V=x5.AD_Data.AD_Avg_V+NoiseFACT*RMS*randn(size(x2.AD_Data.AD_Avg_V));
x6.AD_Data.AD_Avg_V=x6.AD_Data.AD_Avg_V+NoiseFACT*RMS*randn(size(x2.AD_Data.AD_Avg_V));
x7.AD_Data.AD_Avg_V=x7.AD_Data.AD_Avg_V+NoiseFACT*RMS*randn(size(x2.AD_Data.AD_Avg_V));
% plot(x2.AD_Data.AD_Avg_V,'b')
% hold off


write_nel_data('a1002_ABR_1998',x2,0);
write_nel_data('a1003_ABR_1998',x3,0);
write_nel_data('a1004_ABR_1998',x4,0);
write_nel_data('a1005_ABR_1998',x5,0);
write_nel_data('a1006_ABR_1998',x6,0);
write_nel_data('a1007_ABR_1998',x7,0);



