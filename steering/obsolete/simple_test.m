#! /bin/octave -qf

function simple_test(arg = argv())
	
	comid = arg{1}
	
	% Saving the log
	N = datestr(now(), 'yyyymmdd_HHMM');
	dfile = ['.\', comid, '_', N, '_log.txt'];
	diary (dfile);

	
	if size(arg) == 1
		input_dir_arg = "hist";
	else
		input_dir_arg = arg{2};		
	end
	disp(['Input climate location:../', input_dir_arg])

	
	if arg{3}
		LowHighThresh_arg = str2num(arg{3});
	else
		LowHighThresh_arg = 2; % Default value from Fowler was 2
	end
	
	disp(['Threshold for splitting High and Low Freq from IMFs= ', num2str(LowHighThresh_arg)])



end
