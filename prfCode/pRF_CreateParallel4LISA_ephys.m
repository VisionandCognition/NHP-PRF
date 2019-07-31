% CreateParallel Bash scripts for pRF fitting on LISA
m={'lick','aston'};
for i=1:2
    M=m{i};
        joblist.monkey = M;
        joblist.sessions = {...
            'mua', [];...
            'lfp', [];...
            }; % SESSION nWorkers
        joblist.instances = {'1','2','3','4','5','6','7','8'}; % 2 digits, leading zero!
    
    joblist.sessinc     = 1:size(joblist.sessions,1);
    joblist.type        = 'ephys'; % used as label NB! also the folder where data is loaded from
    joblist.modeltype   = 'linear_ephys'; 
    % 'css_hrf' / 'linear_hrf' / 'dog_hrf'
    % 'css_ephys' / 'linear_ephys' / 'dog_ephys'
    
    joblist.xvalmode    = 0; % 0 / 1 / 2
    joblist.resfld      = 'linear_ephys_cv0_neggain';
    
    parallel_fun_dir    = '$TMPDIR/PRF/'; %$TMPDIR is fast 'scratch' space
    parallel_fun        = 'pRF_FitModel_LISA_ephys_neggain';
    
    job_name            = ['FitPRF_' joblist.resfld];
    
    fprintf('\n== Running create_parallel_LISA ==')
    pRF_CreateParallel4LISA_worker_ephys(...
        parallel_fun, joblist, parallel_fun_dir, job_name)
end