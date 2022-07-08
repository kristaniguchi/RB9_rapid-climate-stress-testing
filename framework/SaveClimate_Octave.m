
function SaveClimate_Octave(DataOut, deltaP, deltaT, deltaLowFreqP, deltaSeasonality, info, k, comid)         % later on remove  deltaRRrship_space
    % DataOut: "struct" datatype
	% Timeseries (TS): "table" data type
	% deltaX: double
	
	Out_arg = ['..\out\results\', comid, '\'];
	%mkdir(Out_arg);

	
    if ~info.isOctave
		disp("Use original matlab version of SaveClimateAndFlow.m")
    else
        %disp('Saving .mat binary file.')
		%%% Implemented by Donny Kim. It is a messy approach to save both .mat and .csv files.
		%%% In SCCWRP project, we don't use Wapaba rainfall-runoff model. So, I am getting rid of any Q related components.
		
		
		%% get the data relevant to this deltaRRrship. deltaRRrship = deltaRRrship_space(i)
		%if ~deltaRRrship_space
		%	ThisData = DataOut;
		%else
		%	deltaRRrship = deltaRRrship_space(i); % <- this has to go away later on (DK)
		%	ThisData = DataOut.(['deltaRRrship' num2str(i)]); % <- this has to go away later on (DK)
		%end
		%% DK: Removing this process, as deltaRRrship is no longer used for SCCWRP project.
		%%ThisData = DataOut;

		
		%%% STEP 1: saving .mat binary file
        %%% separate out the items to be saved from the input structure
		ThisData = DataOut;
        TS_P          = table2struct(ThisData.TS_P); 
        TS_T          = table2struct(ThisData.TS_T); 
		
        % remember the perturbations
        perturbations = struct('deltaP', deltaP, 'deltaT', deltaT, 'deltaLowFreqP', deltaLowFreqP, 'deltaSeasonality', deltaSeasonality); 

        % Fowler's original code create name of file using perturbation values.
		% Instead, I will use index (i value) for naming. 
		% Perturbation parameter values will be saved in separate csv file with indexing. This will be done in wrapper loop function.
        OutFile_path = [Out_arg 'StochasticClimate_' ...
            num2str(k) '.mat'];

        % now save: Loop in exmaple_run.m will create csv files, as many as DesiredNumberOfTestRuns (n=1000). Is it desirable?
		if strcmp(fieldnames(ThisData), 'TS_PET') == 1
			TS_PET        = table2struct(ThisData.TS_PET);
			save(OutFile_path, 'TS_P', 'TS_T', 'TS_PET', 'perturbations');
		else
			save(OutFile_path, 'TS_P', 'TS_T', 'perturbations'); %DK: in the future this is the one
		end

		
		
		
		
		%%%% STEP 2: saving Timeseries into .csv file.
		%%%% This part is unavoidably messy, as Octave cannot write neither of Table or Struct into CSV files.
		%%%% Using Cell datatype to write .csv output file.
		%%TS_Q_forModel = preprocess_table2csv(ThisData.TS_Q_ReachInflows, 'TS_Q_forModel');
        %TS_P          = preprocess_table2csv(ThisData.TS_P, 'TS_P'); 
        %TS_T          = preprocess_table2csv(ThisData.TS_T, 'TS_T'); 
        %if ThisData.TS_PET
		%	TS_PET        = preprocess_table2csv(ThisData.TS_PET, 'TS_PET');
		%	TS_merged = [TS_P, TS_T, TS_PET];
		%else
		%	TS_merged = [TS_P, TS_T];
		%end
		%
		%
        %OutFile_path = [Out_arg 'csv\StochasticClimate_' ...
        %    num2str(k) '.csv'];
		%
		%cell2csv([OutFile_path], TS_merged);
		
    end
end




function output = preprocess_table2csv(table, TS_name)
	
	% converting table to cell, and then appending TS_name and VariableNames(colnames) on top two rows of the cell.

	x = table2cell(table);
	y = table.Properties.VariableNames;
	y2 = [repmat({TS_name}, [1, size(y, 2)])];
	
	output = vertcat(y2, y, x);
end