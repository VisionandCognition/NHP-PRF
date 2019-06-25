% CreateParallel Bash scripts for pRF fitting on LISA
%%
joblist.monkey = 'danny';
joblist.sessions = {'medianBOLD_sub-danny', []}; % SESSION nWorkers
%%
%joblist.monkey = 'eddy';
%joblist.sessions = {'medianBOLD_sub-eddy', []}; % SESSION nWorkers


%%
joblist.sessinc = 1:size(joblist.sessions,1); 
%joblist.slicechunks = {'01:14','15:28','29:42','43:56'}; % 2 digits, leading zero!
joblist.slicechunks = {'01:15','16:30','31:45','46:60'}; % 2 digits, leading zero!
joblist.type = 'avg';
joblist.hrf = ''; %'defaultHRF';

parallel_fun_dir    = '$TMPDIR/PRF/'; %$TMPDIR is fast 'scratch' space
parallel_fun        = 'pRF_FitModel_LISA_avg';
job_name            = 'FitPRF_avg';
 
disp('== Running create_parallel_LISA ==')
 
pRF_CreateParallel4LISA_worker_avg(...
    parallel_fun, joblist, parallel_fun_dir, job_name)
 