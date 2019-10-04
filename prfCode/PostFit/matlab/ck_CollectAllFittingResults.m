%% MRI ====================================================================
monkeys = {'eddy','danny'};

models = {...
    'linhrf_cv1_dhrf','linhrf_cv1_mhrf',...
    'linhrf_cv1_dhrf_neggain','linhrf_cv1_mhrf_neggain',...
    'doghrf_cv1_dhrf','doghrf_cv1_mhrf',...
    'csshrf_cv1_dhrf','csshrf_cv1_mhrf'};

ck_GetMRI_pRF(monkeys,models);


%% EPHYS ==================================================================
monkeys = {'lick','aston'};

models = {...
    'linear_ephys_cv1',...
    'linear_ephys_cv1_neggain',...
    'dog_ephys_cv1',...
    'css_ephys_cv1',...
    'classicRF'};

ck_GetEphys_pRF(monkeys,models);


%% COMBINE ================================================================
fitres_path = ...
    '/Users/chris/Dropbox/CURRENT_PROJECTS/NHP_MRI/Projects/pRF/FitResults/';

% load mri
fprintf('Loading MRI\n');
R_MRI = load(fullfile(fitres_path,'MRI','Combined','AllFits_MRI_cv1'));
% load ephys
fprintf('Loading EPHYS\n');
R_EPHYS = load(fullfile(fitres_path,'ephys','Combined','AllFits_ephys_cv1'));

% save
fprintf('Saving results from both modalities\n');
save(fullfile(fitres_path,'MultiModal','AllFits_cv1'),'R_MRI','R_EPHYS');

%% PERFORM SOME BASIC PROCESSING ON COLLECTED RESULTS =====================



%% CREATE LARGE TABLES WITH ALL RESULTS ===================================