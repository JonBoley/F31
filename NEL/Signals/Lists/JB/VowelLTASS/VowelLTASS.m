% Written by GE_MH 04nov2003. Modified EDY 5/10/05.
% Modified M Heinz 8/3/07
% Modified M Heinz 12/15/08

global signals_dir

files ={
'BBN_A_.wav',
'BBN_A_N.wav',
'BoyFell_.wav',
'BoyFell_N.wav'};

Lchannel.file_list = cell(length(files),1);

ii=0;
for i=1:length(files)
    ii=ii+1;
    Lchannel.file_list{ii} = sprintf('%sMH\\ISH2009\\%s', signals_dir,files{i});
end
