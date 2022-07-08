#! /bin/octave -qf
% Modified by Donny Kim. Last update: June, 2022
% Based on January 2021 by Keirnan Fowler, University of Melbourne, fowler.k@unimelb.edu.au

% Integrated framework for rapid climate stress testing on a monthly timestep
% by Keirnan Fowler, Natasha Ballis, Avril Horne, Andrew John, Rory Nathan and Murray Peel
% Licence: CC BY 3.0 - see https://creativecommons.org/licenses/by/3.0/au/



% Following function is designed to be called from R.
% R will pass A SINGLE argument (e.g. "comid_20332660") to run the function.
% Paralell processing is purely governed by R using "parallel", "doParallel", and "foreach" package. This Octave code is not implemented with any parallelization by itself.


%if size(arg) == 1
%	comid_arg = arg{end};
%	disp('Default input location: ../hist')
%else
%	comid_arg = arg{1};
%	input_dir_arg = arg{end};
%	disp(['Input location argument detected: ../', input_dir_arg])
%end


comid = 'comid_10000042'; % Just in case
input_dir_arg = 'hist_raw';
mkdir('..\out\');

% add path to directories with framework code
% Equivalent to 'sourcing' functions in R
addpath('..\framework\');
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0; % check if we are in Octave
if isOctave, addpath(genpath('..\framework\octave\')); pkg load statistics; %pkg load tablicious;
end


% importing climatic hist data.
% * Weirdly, "readtable" is custom function made by Fowler et al., and it cannot read csv file if it has a column name that only consist of numbers.
% * This is the reason why using input argument such as "comid_20332660" instead of simply doing "20332660".
precip_monthly = readtable(['..\', input_dir_arg, '\precip_monthly_', comid(7:end), '.csv']); % comid(7:end) extracts actual comid (e.g. 20332660) from into argument string such as "comid_20332660"
T_monthly = readtable(['..\', input_dir_arg, '\T_monthly_', comid(7:end), '.csv']);
%PET_monthly = readtable('hist\PET_monthly.csv'); % DK: modified framework is capable of handling PET input. It is optional input at this time.


% DK: flow data will not be used for SCCWRP project. Any flow related components are completely removed below this section.
%flow_monthly.RepresentativeCatchments_mm = readtable('hist\flow_monthly_representative_catchments.csv');
%flow_monthly.ReachInflows_ML = readtable('hist\flow_monthly_reach_inflows.csv');


% Store P and T in a single structure
%HistoricData = struct('precip_monthly', precip_monthly, 'T_monthly', T_monthly, 'PET_monthly', PET_monthly);
HistoricData = struct('precip_monthly', precip_monthly, 'T_monthly', T_monthly);


% build a info file, using simple function. Note that any of "flow" related components are removed.
info = comid_info(comid);  % specify information and settings, including the following:
						% info.SubareaDetails         Subarea names, areas, and weightings to use in aggregation of Intrinsic Mode Functions (IMFs) (see function AggregateIMFs)
						% info.pars                   Other general settings relevant to stochastic climate generation and the way the data is used
info.isOctave = isOctave;


% aggregate all data to annual
% note, the results are appended to the existing data structure 'HistoricData'
HistoricData = AggregateToAnnual(HistoricData, info.pars);


% Removing variables other than HistoricData and info for memory management... but it is unnecessary.
clearvars -except HistoricData info comid


%% Split precipitation into high and low frequency using Empirical Mode Decomposition
% use the CEEMDAN algorithm (Complete Ensemble Empirical Mode Decomposition with Adaptive Noise)
% note, the results are appended to the existing data structure 'HistoricData'
disp('Finished loading data.  Starting empirical mode decomposition using CEEMDAN (it takes a WHILE)')

% !!!!!!!!
% !!! This is the SINGLE USER DEFINED PARAMETER (not really, but technically)
% set the threshold that determines which IMFs are 'low' and 'high' (see paper, Section ##)
info.LowHighThresh = 1; % IMFs 1 and 2 are high frequency, everything else is low
% !!!!!!!!

% run CEEMDAN
% DK: Currently getting all the annual low_freq as 0. This may be due to parameter setting OR short historic data timeperiod
HistoricData = split_high_low_using_CEEMDAN(HistoricData, info); % this step may randomly fail. I am getting LowFreq as 0 all the time.

%% Conduct pre-analysis to inform stochastic generation of the low-frequency component
disp('Done.  Starting pre-analysis for low-frequency component.');
info.LowFreq_PreAnalysis_Outputs = LowFreq_PreAnalysis(HistoricData, info);

%%% Saving HistoricData and info, before perturbation
%%% Note that this
%addpath(genpath('..\framework\octave\')) % we already did it, but just to make sure.
save_IntermediateFiles_Clim(HistoricData, info, comid);%%% Loading saved intermediate .mat binary files

% Loading function for intermediate outputs.
[HistoricData, info] = load_IntermediateFiles_Clim(comid);





%%% testing purpose only =================================================================================================
%samples = info.SubareaList';
%temp.A = struct2table(info.LowFreq_PreAnalysis_Outputs.A);
%for i = 1:size(samples, 2);
%	temp.(samples{i}) = struct((info.LowFreq_PreAnalysis_Outputs.(samples{i})));
%  temp.(samples{i}) = struct2table(temp.(samples{i}));
%end
%temp.TS_rand = info.LowFreq_PreAnalysis_Outputs.TS_rand;
%%%%%% testing endes here===============================================================================================

%% case 2: full stochastic dataset
%  generate a full dataset of stochastic data covering the entire
%  stress testing space, as defined by the axis limits and
%  gradations specified in Table # in the paper.  Whereas the code above
%  placed all the outputs into a single structure (StochPertData), here
%  there are too many so we save them as individual .mat files.

% specify axis gradations
%deltaP_space           = [-0.40 -0.35 -0.30 -0.25 -0.20 -0.15 -0.10 -0.05 0 .05 .10 .15]; % rainfall proportional change "increase 0.1
deltaP_space           = [-0.10 -0.05 0 .05 .10];
%deltaT_space           = [0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0];                             % additional degrees of warming 2.0
deltaT_space           = [0 0.5 1.0 1.5 2.0];

%deltaLowFreqP_space    = [-0.03 -0.015 0 +0.015 +0.03 +0.045 +0.06 +0.075];               % changes n Hurst Coefficient %0
deltaLowFreqP_space    = [0];               % changes n Hurst Coefficient %0
%deltaSeasonality_space = [-0.06 -0.03 0.00 +0.03 +0.06 +0.09 +0.12 +0.15];                % changes in seasonality %0
deltaSeasonality_space = [0];                % changes in seasonality %0
%deltaRRrship_space     = [-50 -43.75 -37.5 -31.25 -25 -18.75 -12.5 -6.25 0 6.25 12.5];    % shift in rainfall runoff relationship
%deltaRRrship_space = [];



%%% Main stochastic time-series generation & perturbation loop
k = 0; % index to match perturbation parameter-set and TS_P & TS_T later on
perturbations_save = {};
TS_P_save = {};
TS_T_save = {};
mkdir(['..\out\results\', comid, '\']);


for deltaP = deltaP_space
	for deltaT = deltaT_space
		for deltaLowFreqP = deltaLowFreqP_space
			for deltaSeasonality = deltaSeasonality_space

				DesiredNumberOfTestRuns = 25; % Setting this number to match with perturbation space numbers to get full results.
				if rand() < (DesiredNumberOfTestRuns / (5*5*1*1)) % DK: this controls the random sampling. What an odd way to do so...

					% Update index k
					k = k+1;

					% run the stochastic data routines
					disp([comid, " StochPert: ", num2str(k), " out of ", num2str(DesiredNumberOfTestRuns)]);
					DataOut = GetStochPertData(deltaP, deltaT, deltaLowFreqP, deltaSeasonality, HistoricData, info); %DK: using version of GetStochPertData where any flow/rainfall-runoff related components are disabled.

					% DK: Recording TS_P, TS_T, and perturbation parameters for given index k.
					temp = {k, deltaP, deltaT, deltaLowFreqP, deltaSeasonality};
					perturbations_save = vertcat (perturbations_save, temp);
					TS_P_save = [TS_P_save, table2cell(DataOut.TS_P)(:,3)];
					TS_T_save = [TS_T_save, table2cell(DataOut.TS_T)(:,3)];

					%SaveClimate_Octave(DataOut, deltaP, deltaT, deltaLowFreqP, deltaSeasonality, info, k, comid); %Optional when you want to save .mat binary file. DK: I personally think there's no point of saving it. Waste of disk space.

				end
			end
		end
	end
end


% Retouching recorded TS_P and TS_T
colnames = {'Year', 'Month'};
for i = 1:k
	colnames = [colnames, i];
end

TS_P_save = [table2cell(DataOut.TS_P)(:,1:2), TS_P_save]; TS_P_save = vertcat(colnames, TS_P_save);
TS_T_save = [table2cell(DataOut.TS_T)(:,1:2), TS_T_save]; TS_T_save = vertcat(colnames, TS_T_save);

% Save to CSV files
cell2csv (['..\out\results\', comid ,'\TS_P.csv'], TS_P_save);
cell2csv (['..\out\results\', comid ,'\TS_T.csv'], TS_T_save);
cell2csv (['..\out\results\', comid ,'\perturbations_save.csv'], perturbations_save);


% To extract Low_Freq_Comp for non-perturbed stochastically generated TS. Inteded for some sort of visual inspection similar to Fig 4 (Fowler et al., 2022).
[~, extra] = GetStochPertData(0, 0, 0, 0, HistoricData, info);
TS_P_ann_LowFreq = table2cell(extra.extra_b.TS_P_ann_LowFreq);
TS_P_ann_HighFreq = table2cell(extra.TS_P_ann_HighFreq);
TS_P_ann_HiLoFreq = [num2cell(1:info.pars.StochRepLen_yrs)', num2cell(extra.extra_f.years)', TS_P_ann_LowFreq, TS_P_ann_HighFreq];
TS_P_ann_HiLoFreq = vertcat({'Year', 'RandomizedFrom', 'LowFreq', 'HighFreq'}, TS_P_ann_HiLoFreq);

cell2csv (['..\out\results\', comid ,'\TS_HiLowFreq.csv'], TS_P_ann_HiLoFreq);


%csvwrite(['..\out\intermediate\', comid ,'\IMFs.csv'], HistoricData.precip_annual_IMFs.comid_10000042)
