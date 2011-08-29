function plot_dps

global dp_data_dir han z animal date

plot(han.dp_curve,z.DpoaeData(:,3),z.DpoaeData(:,4),'*k:',z.Dpoaefreqs(1,:),z.DpoaeSpectra(15,:),'r:')
	
	