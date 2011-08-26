function abr_setup

host = lower(getenv('hostname'));
switch (host)
case {'south-chamber'}
	addpath c:\Users\GE\Matlab_ABR;
   addpath c:\Users\GE\Matlab_ABR\file_manager;
   addpath c:\Users\GE\Matlab_ABR\signal_averager;
end   
