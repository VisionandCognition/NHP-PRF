CV='cv0'; %'cv0'/'cv1'

Do.MRI = false;
Do.EPHYS = true;
Do.COMBINE = true;

%% MRI ====================================================================
if Do.MRI
    monkeys = {'eddy','danny'};
    
    models = {...
        ['linhrf_' CV '_dhrf'],['linhrf_' CV '_mhrf'],...
        ['linhrf_' CV '_dhrf_neggain'],['linhrf_' CV '_mhrf_neggain'],...
        ['doghrf_' CV '_dhrf'],['doghrf_' CV '_mhrf'],...
        ['csshrf_' CV '_dhrf'],['csshrf_' CV '_mhrf']};
    
    output = ['AllFits_MRI_' CV];
    
    ck_GetMRI_pRF(monkeys,models,output);
end

%% EPHYS ==================================================================
if Do.EPHYS
    monkeys = {'lick','aston'};
    
    models = {...
        ['linear_ephys_' CV],...
        ['linear_ephys_' CV '_neggain'],...
        ['dog_ephys_' CV],...
        ['css_ephys_' CV],...
        'classicRF'};
    output = ['AllFits_ephys_' CV];
    
    ck_GetEphys_pRF(monkeys,models,output);
end

%% COMBINE ================================================================
if Do.COMBINE
    fitres_path = ...
        '/Users/chris/Dropbox/CURRENT_PROJECTS/NHP_MRI/Projects/pRF/FitResults/';
    
    % load mri
    fprintf('Loading MRI\n');
    R_MRI = load(fullfile(fitres_path,'MRI','Combined',['AllFits_MRI_' CV]));
    R_MRI = R_MRI.R;
    % load ephys
    fprintf('Loading EPHYS\n');
    R_EPHYS = load(fullfile(fitres_path,'ephys','Combined',['AllFits_ephys_' CV]));
    R_EPHYS = R_EPHYS.R;
    
    % save
    fprintf('Saving results from both modalities\n');
    save(fullfile(fitres_path,'MultiModal',['AllFits_' CV]),'R_MRI','R_EPHYS');
end