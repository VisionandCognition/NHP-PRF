% CreateParallel Bash scripts for pRF fitting on LISA
m={'danny','eddy'};
for i=1:2
    M=m{i};
    if strcmp(M,'danny')
        joblist.monkey = 'danny';
        joblist.sessions = {...
            'AllSessions-avg-cv0', [];...
            }; % SESSION nWorkers
        joblist.slicechunks = {'01:14','15:28','29:42','43:56'}; % 2 digits, leading zero!
    elseif strcmp(M,'eddy')
        joblist.monkey = 'eddy';
        joblist.sessions = {...
            'AllSessions-avg-cv0', [];...
            }; % SESSION nWorkers
        joblist.slicechunks = {'01:15','16:30','31:45','46:60'}; % 2 digits, leading zero!
    end
    
    joblist.sessinc     = 1:size(joblist.sessions,1);
    joblist.type        = 'cv0'; % used as label NB! also the folder where data is loaded from
    joblist.hrf         = 'defaultHRF'; % 'HRF_monkey' / 'defaultHRF' / 'none'
    joblist.modeltype   = 'linear_hrf'; 
    % 'css_hrf' / 'linear_hrf' / 'dog_hrf'
    % 'css_ephys' / 'linear_ephys' / 'dog_ephys'
    
    joblist.xvalmode    = 0; % 0 / 1 / 2
    joblist.resfld      = 'linhrf_cv0_dhrf';
    
    parallel_fun_dir    = '$TMPDIR/PRF/'; %$TMPDIR is fast 'scratch' space
    parallel_fun        = 'pRF_FitModel_LISA_cv';
    
    job_name            = ['FitPRF_' joblist.resfld];
    
    fprintf('\n== Running create_parallel_LISA ==')
    pRF_CreateParallel4LISA_worker_cv0(...
        parallel_fun, joblist, parallel_fun_dir, job_name)
end