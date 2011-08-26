function update_textfile_format(fname)
%

% AF 12/11/01

global spikes

[pa,ff] = fileparts(fname);
if (exist(ff)==2)
   x = eval(ff,['Can''t evaluate ''' ff '''']);
else
  return
end
if (isfield(x,'Stimuli') & isfield(x,'Line') & isfield(x,'General') & isfield(x,'spikes'))
   if (~isfield(x.General,'trigger'))
      x.General.Trigger = '';
   end
   if (~isfield(x.General,'run_errors'))
      x.General.run_errors = '';
   end
   if (~isfield(x.Stimuli,'fully_presented_lines'))
      if (~isempty(x.Line))
         tmp = struct2cell(x.Line.attens);
         x.Stimuli.fully_presented_lines = size(tmp{1},1);
      else
         x.Stimuli.fully_presented_lines = 0;
      end
   end
   if (~isfield(x.Stimuli,'fully_presented_stimuli'))
      if (~isempty(x.Line))
         tmp = struct2cell(x.Line.attens);
         x.Stimuli.fully_presented_stimuli = size(tmp{1},1);
      else
         x.Stimuli.fully_presented_stimuli = 0;
      end
   end
   %% Filling the global 'spikes' variable.
   nChannels = length(x.spikes); 
   spikes.last  = zeros(1,nChannels);
   spikes.times = x.spikes;
   for ii = 1:nChannels
      spikes.last(ii) = size(spikes.times{ii},1);
   end
   %
   rc = write_nel_data(fname,x,1);
   if (rc < 0)
      error(['problem in saving file ''' fname '''']);
   end
end

   