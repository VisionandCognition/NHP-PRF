% CreateParallel Bash scripts for pRF fitting on LISA
m={'danny','eddy'};
for i=1:2
    M=m{i};
    if strcmp(M,'danny')
        joblist.monkey = 'danny';
        joblist.sessions = {...
            'AllSessions-only_avg', [];...
            'AllSessions-avg-even', [];...
            'AllSessions-avg-odd', [];...
            }; % SESSION nWorkers
        joblist.slicechunks = {'01:14','15:28','29:42','43:56'}; % 2 digits, leading zero!
    elseif strcmp(M,'eddy')
        joblist.monkey = 'eddy';
        joblist.sessions = {...
            'AllSessions-only_avg', [];...
            'AllSessions-avg-even', [];...
            'AllSessions-avg-odd', [];...
            }; % SESSION nWorkers
        joblist.slicechunks = {'01:15','16:30','31:45','46:60'}; % 2 digits, leading zero!
    end
    
    joblist.sessinc = 1:size(joblist.sessions,1);
    joblist.type = 'avg';
    joblist.hrf = 'HRF_monkey'; %'defaultHRF';
    
    parallel_fun_dir    = '$TMPDIR/PRF/'; %$TMPDIR is fast 'scratch' space
    parallel_fun        = 'pRF_FitModel_LISA_avg';
    job_name            = 'FitPRF_avg';
    
    fprintf('\n== Running create_parallel_LISA ==')
    
    pRF_CreateParallel4LISA_worker_avg(...
        parallel_fun, joblist, parallel_fun_dir, job_name)
end