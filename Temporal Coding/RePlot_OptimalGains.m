%% Test OptimalGain
impairment = 'mixed';
phone = 12;
levels = 60:10:100;
gains = 40:-5:-40;
note = 'Apr_06_09';

[gain1,gain2,err1,err2] = OptimalGain2(levels,gains,sprintf('archive\\phone%d\\',phone),impairment,note);

