% ck_MakePlots
% Takes the pre-processed fitting result tables and creates comparison plots

%% Paths ==================================================================
BaseFld = pwd;
ResFld = ...
    '/Users/chris/Documents/MRI_ANALYSIS/NHP-PRF/FitResults/MultiModal';
T='Tables_max';

%% Load ===================================================================
fprintf('Loading results table. Please be patient, this will take a while..\n');
tic; load(fullfile(ResFld,T)); t_load=toc;
fprintf(['Loading took ' num2str(t_load) ' s\n']);

