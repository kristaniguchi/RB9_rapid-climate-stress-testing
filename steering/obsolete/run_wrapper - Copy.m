% Created by Donny Kim. Last update: June 26th, 2022
% Created January 2021 by Keirnan Fowler, University of Melbourne, fowler.k@unimelb.edu.au

% Integrated framework for rapid climate stress testing on a monthly timestep
% by Keirnan Fowler, Natasha Ballis, Avril Horne, Andrew John, Rory Nathan and Murray Peel

% Licence: CC BY 3.0 - see https://creativecommons.org/licenses/by/3.0/au/

% Use this file to run the example given in the paper.  Specifically, this
% file can be used to demonstrate for the following two cases:
% - case 1: only generate the stochastic data required to create the
%           plots in the paper, and no more.  Then do the plotting.
% - case 2: generate a full dataset of stochastic data covering the entire
%           stress testing space, as defined by the axis limits and
%           gradations specified in Table # in the paper

% This file uses the code contained in folder '\framework', producing
% stochastic output saved in folder '\example\out' (in .mat format), then
% using code from example_plotting.m to create plots.

clear; close all;

% add path to directories with framework code
addpath('..\framework\');
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0; % check if we are in Octave
if isOctave, addpath(genpath('..\framework\octave\')); pkg load statistics; pkg load tablicious;
end
% start timing the code
tic;

%% load inputs including historic data

% climatic data
addpath(genpath('..\framework\octave\'));
precip_monthly = readtable('..\hist\precip_monthly_20325361.csv');
T_monthly = readtable('..\hist\T_monthly_20325361.csv');
comid = 20325361
%PET_monthly = readtable('hist\PET_monthly.csv');



% flow data - note, this is not used directly by the framework but is used
% to produce plots that appear in the paper.  The unit of flow varies with
% context.
%flow_monthly.RepresentativeCatchments_mm = readtable('hist\flow_monthly_representative_catchments.csv');
%flow_monthly.ReachInflows_ML = readtable('hist\flow_monthly_reach_inflows.csv');

% store all this data in a single structure
%HistoricData = struct('precip_monthly', precip_monthly, 'T_monthly', T_monthly, 'PET_monthly', PET_monthly, 'flow_monthly', flow_monthly);
HistoricData = struct('precip_monthly', precip_monthly, 'T_monthly', T_monthly);


info = comid_info('comid_20325361');  % specify information and settings, including the following:
                        % info.WapabaParSets          Parameter sets for pre-calibrated WAPABA rainfall runoff model - see Supplementary Material ##
                        % info.SubareaDetails         Subarea names, areas, and weightings to use in aggregation of Intrinsic Mode Functions (IMFs) (see function AggregateIMFs)
                        % info.RepCatchDetails        Each subarea has a "representative" catchment.  This table stores the characteristics of these catchments.  See paper Section ## and Supplementary Material ##.
                        % info.FlowConversionFactors  Each reach of the Goulburn River receives inflows from multiple subareas.  This table of factors describes how.  See paper Section ## and Supplementary Material ##.
                        % info.pars                   Other general settings relevant to stochastic climate generation and the way the data is used

% aggregate all data to annual
% note, the results are appended to the existing data structure 'HistoricData'
HistoricData = AggregateToAnnual(HistoricData, info.pars);

info.isOctave = isOctave;
%clearvars -except HistoricData info

%% Split precipitation into high and low frequency using Empirical Mode Decomposition
% use the CEEMDAN algorithm (Complete Ensemble Empirical Mode Decomposition with Adaptive Noise)
% note, the results are appended to the existing data structure 'HistoricData'
disp('Finished loading data.  Starting empirical mode decomposition using CEEMDAN...')
disp('Running CEEMDAN algorithm takes A LOT of time. So be patient and wait for it...')

% set the threshold that determines which IMFs are 'low' and 'high' (see paper, Section ##)
info.LowHighThresh = 2; % IMFs 1 and 2 are high frequency, everything else is low

% run CEEMDAN
HistoricData = split_high_low_using_CEEMDAN(HistoricData, info);

%% Conduct pre-analysis to inform stochastic generation of the low-frequency component
disp('Done.  Starting pre-analysis for low-frequency component.');
info.LowFreq_PreAnalysis_Outputs = LowFreq_PreAnalysis(HistoricData, info);


%%% Saving HistoricData and info, before perturbation
%%% Note that this
addpath(genpath('..\framework\octave\')) % we already did it, but just to make sure.
%save_IntermediateFiles_FULL(HistoricData, info) % This function is included in /framework/octave directory


%%% Loading saved intermediate .mat binary files
[HistoricData, info] = load_IntermediateFiles_FULL('out\intermediate\');

timing.init = toc;



%%% testing purpose only =================================================================================================
%samples = info.SubareaList';
%temp.A = struct2table(info.LowFreq_PreAnalysis_Outputs.A);
%for i = 1:size(samples, 2);
%	temp.(samples{i}) = struct((info.LowFreq_PreAnalysis_Outputs.(samples{i})));
%  temp.(samples{i}) = struct2table(temp.(samples{i}));
%end
%temp.TS_rand = info.LowFreq_PreAnalysis_Outputs.TS_rand;
%%%%%% testing endes here===============================================================================================

tic;
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
%deltaRRrship_space     = [-15 0];    % shift in rainfall runoff relationship

% generate data
addpath('..\framework\');
i = 1;
perturbations_save = {};

% test==========================================================================
i = 1;
perturbations_save = {};

DataOut = GetStochPertData(deltaP_space(1), deltaT_space(1),...
    deltaLowFreqP_space(1), deltaSeasonality_space(1), deltaRRrship_space, HistoricData, info); %DonnyKim: May be we could try to get rid of deltaRRrship_space afterall.).

prettyprint(DataOut.TS_P (1:15, :));

deltaRRrship = deltaRRrship_space(i);
if ~deltaRRrship
	ThisData = DataOut;
else
	deltaRRrship = deltaRRrship_space(i); % <- this has to go away later on (DK)
	ThisData = DataOut.(['deltaRRrship' num2str(1)]); % <- this has to go away later on (DK)
end


TS_P          = preprocess_table2csv(ThisData.TS_P, 'TS_P');
TS_T          = preprocess_table2csv(ThisData.TS_T, 'TS_T');
TS_PET        = preprocess_table2csv(ThisData.TS_PET, 'TS_PET');

TS_merged = [TS_P, TS_T];

OutFile_path = ['..\out\out_data\csv\' 'StochasticClimate_' num2str(i) '.csv'];

cell2csv([OutFile_path], TS_merged);

% end test =====================================================================



% test2 ========================================================================
addpath('..\framework\');

[DataOut, extra] = GetStochPertData(deltaP_space(1), deltaT_space(1), deltaLowFreqP_space(1), deltaSeasonality_space(1), HistoricData, info);
% end test =====================================================================

addpath('..\framework\');
k = 1; % index
perturbations_save = {};
OutFile_path = ['..\out\out_data\];
mkdir('..\out\out_data\csv');

for deltaP = deltaP_space
    for deltaT = deltaT_space
        for deltaLowFreqP = deltaLowFreqP_space
            for deltaSeasonality = deltaSeasonality_space

                DesiredNumberOfTestRuns = 10;
                if rand() < (DesiredNumberOfTestRuns / (5*5)) % DK: this controls the random sample

                    % run the stochastic data routines
                    disp([num2str(k), " out of target number: ", num2str(DesiredNumberOfTestRuns)]);
                    DataOut = GetStochPertData(deltaP, deltaT, deltaLowFreqP, deltaSeasonality, deltaRRrship_space, HistoricData, info); %DonnyKim: May be we could try to get rid of deltaRRrship_space afterall.).
                    disp('Perturbation complete, saving .mat binary and CSV files');

					          % DK: Saving selected perturbation parameters with indexing. It will be a single CSV file.
					          temp = {k, deltaP, deltaT, deltaLowFreqP, deltaSeasonality};
					          perturbations_save = vertcat (perturbations_save, temp);

                    for i = 1:size(deltaRRrship_space, 2)
                        SaveClimateAndFlow(i, DataOut, deltaP, deltaT, deltaLowFreqP, deltaSeasonality, deltaRRrship_space, info, k, OutFile_path, comid);
                    end
                    % save to file. DK: in our case, both .mat and .csv files (total number of files = i*2)
                    k = k+1;
                    % *note, to save disc space the user could alter this to only save one copy of the climate inputs at this point, and then
                    %  a separate file for each deltaRRrship.

                end
            end
        end
    end
end
cell2csv ("perturbations_save.csv", perturbations_save); % Writing perturbation parameters into actual CSV.

timing.case2_stochgen = toc;

%% report back on timing
disp('All done.  Run times were:')
disp('Initialisation including data loading and low frequency pre-analysis: ')
disp([num2str(timing.init) ' seconds. ']);
disp('Stochastic generation for case 1')
disp([num2str(timing.case1_stochgen) ' seconds. ']);
disp('Plotting for case 1')
disp([num2str(timing.plotting) ' seconds. ']);
disp('Stochastic generation for case 2')
disp([num2str(timing.case2_stochgen) ' seconds. ']);
