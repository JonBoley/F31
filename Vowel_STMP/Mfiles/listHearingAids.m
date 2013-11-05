%%
dirList = dir('*WAVreBFi*');
for ii=1:numel(dirList)
    fprintf('%s',dirList(ii).name(1:end-2));
    eval(sprintf('x=%s;',dirList(ii).name(1:end-2)));
    if strncmpi(x.Stimuli.Condition.HearingAid_Nonlinear,'y',1)
        fprintf(' Nonlinear');
    end
    if strncmpi(x.Stimuli.Condition.HearingAid_Linear,'y',1)
        fprintf(' Linear');
    end
    if strncmpi(x.Stimuli.Condition.HearingAid_None,'y',1)
        fprintf(' None');
    end
    fprintf('\n');
end
