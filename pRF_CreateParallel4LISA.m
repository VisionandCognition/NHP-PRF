% CreateParallel Bash scripts for pRF fitting on LISA
joblist.monkey = 'danny';
joblist.sessions = {...
    '20171116',... 1
    '20171129',... 2
    '20171207',... 3
    '20171214',... 4
    '20171220',... 5
    '20180117',... 6
    '20180124',... 7
    '20180125',... 8
    '20180131',... 9
    '20180201'...  10
    };%[1:10]; % SESSION
joblist.sessinc = 1:10; 

parallel_fun_dir = '"$TMPDIR"/PRF/'; %$TMPDIR is fast 'scratch' space
parallel_fun = 'parallelFun_LISA';
job_name = 'FitPRF_PerSession';

disp('== Running create_parallel_LISA ==')
create_parallel_LISA(parallel_fun, joblist, parallel_fun_dir, job_name)
