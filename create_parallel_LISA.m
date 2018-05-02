function create_parallel_LISA(parallel_fun, joblist, parallel_fun_dir, job_name)
% This function serves to create parallel Jobs for a given script
% to parallalise it for a given list of jobs.
%
% Therefore, it creates executable sh-scripts, which execute each single
% job, as well as one script which passes all the single job scripts to
% the cluster computer. (LISA @ SurfSara).
%
% This version requsts 600MB free on the host machine. If your job does not
% require free memory (also matlab already needs quite a bit), use
% create_parallel_nomem.m.
%
% IMPORTANT: The function runs on linux: all pathes must be linux style
%
% PARAMETERS
%   parallel_fun: name of the function (e.g. 'example_function', NOT
%               'example_script.m', 'example_function()', etc.). This
%               function must take exactly 1 argument as input. For an
%               example, see parallel_example_fun.m;
%
%   joblist: a structure with information on what jobs to create
%               [this was a 1xn vector of elements which will be one after
%                the other passed to the function parallel_fun.]
%
%   parallel_fun_dir: path to parallel_fun, will be used to add to the
%               matlab path on the remote machines (the nathans)
%
% OPTIONAL
%   job_name: string to identify your job. Subdirectories with this name will
%       be created inside the specified batch directory and the log
%       directory, so that multiple jobs can be executed at the same time.
%       If no jobname is given, the name of prallel_fun together with the
%       current date and time will be used as job name.

%% basic input checks -----------------------------------------------------
% check if parallel_fun ends with .m or ()
if length(parallel_fun) >= 2 && (strcmp(parallel_fun(end-1:end), '.m') || ...
        strcmp(parallel_fun(end-1:end), '()'))
    parallel_fun = parallel_fun(1:end-2);
end

%% set parameters ---------------------------------------------------------
% default name for jobs
if ~exist('job_name', 'var')
    job_name = [parallel_fun '_' datestr(now, 'yyyymmddTHHMMSS')];
end

% set directory where the sh-files should be created
% will be created, if it does not exist
project_dir = '/home/pcklink/PRF'; % must be the ABSOLUT path

% batch files will be written to job_name/batchdirs
% batch_dir = '~/condor/batch_dir/'  % must be the ABSOLUT path
batch_dir = [project_dir '/' job_name '/batch_dir']; % add jobname

% set log folder
% will be created, if it does not exist
log_file_dir = [project_dir '/' job_name '/logs/']; % add jobname

%% location of scripts ----------------------------------------------------
% set location of execute_matlab_process.sh
execute_matlab_process_sh = '"$TMPDIR"/BashScripts/execute_Compiledmatlab_process_SARA.sh'; % must be the ABSOLUTE path
% set location of execute_matlab_process_sh
generate_submit = '"$TMDPIR"/BashScripts/SubmitterOfAll.sh'; % must be the ABSOLUT path

%% PROCESSING STARTS FROM HERE (no more parameters to check) ==============
%% create batch & log folder ----------------------------------------------
disp('Creating batch & log folders')
[success, message] = mkdir(batch_dir);
if ~success
    error(['Could not create directory for batch_dir: ' message])
end

if ispc
    error('Windows will not work due to paths. Please run from a Linux machine')
else
    [success, message] = mkdir(log_file_dir);
    if ~success
        error(['Could not create directory for log_file_dir: ' message])
    end
end

%% Check if main batch files already exists -------------------------------
% set initial behaviour, if you want to overwrite sh-files which are there
overwrite_file = 'ask';
disp('Creating batch files')

% check if main sh-file to start all jobs exists
filename_all = sprintf('send_all_bootstrap_jobs.sh');
fullfilename_all = [batch_dir '/' filename_all];
if exist(fullfilename_all, 'file')
    disp(' ')
    disp(['File ' fullfilename_all ' already exist.'])
    overwrite_file = input('Should it be overwritten? [y, n, a (all)]: ', 's');
    if ~(strcmpi(overwrite_file, 'y') || strcmpi(overwrite_file, 'a'))
        error(['File ' filename_all ' already exists and should not be overwritten. Solve problem and start again.'])
    end
    delete(fullfilename_all)
end

filename_ori = sprintf('send_all_ori_jobs.sh');
fullfilename_ori = [batch_dir '/' filename_ori];
if exist(fullfilename_ori, 'file')
    disp(' ')
    disp(['File ' fullfilename_ori ' already exist.'])
    overwrite_file = input('Should it be overwritten? [y, n, a (all)]: ', 's');
    if ~(strcmpi(overwrite_file, 'y') || strcmpi(overwrite_file, 'a'))
        error(['File ' filename_ori ' already exists and should not be overwritten. Solve problem and start again.'])
    end
    delete(fullfilename_ori)
end

%% Create the batch files -------------------------------------------------
% The main batch file handles passing the single job batch files to the
% cluster (using lisa_generate_submit.sh)
display(['Creating main batch file: ' fullfilename_all])
fid_commit_all = fopen(fullfilename_all, 'w');

% ensure that the right shell is used !#/bin/sh
fprintf(fid_commit_all, '#!/bin/bash\n');
% add comment that THIS file submits the stuff to condor
fprintf(fid_commit_all, '#\n');
fprintf(fid_commit_all, '# This bash-script submits all jobs to the server, instead of running them locally.\n');
fprintf(fid_commit_all, ['# If you want to submit only some jobs to the server, simply add a "#" in front of \n' ...
    '#the ones you like to ommit and execute the script then.\n']);
fprintf(fid_commit_all, '#\n');

% [CK] Don't think I need this... ======
% fid_commit_ori = fopen(fullfilename_ori, 'w');
% fprintf(fid_commit_ori, '#!/bin/bash\n');
% % add comment that THIS file submits the stuff to condor
% fprintf(fid_commit_ori, '#\n');
% fprintf(fid_commit_ori, '# This bash-script submits all jobs to the server, instead of running them locally.\n');
% fprintf(fid_commit_ori, ['# If you want to submit only some jobs to the server, simply add a "#" in front of \n' ...
%     '# the ones you like to ommit and execute the script then.\n']);
% fprintf(fid_commit_ori, '#\n');
% ======================================

% create all single job batchfiles, and add for each a call in the main
% batch file
for job_ind = 1:length(joblist.sessions)
    %% create batchfile for current job -----------------------------------
    % create/overwrite file
    filename = sprintf('run_job_Ses-%s_%s.sh', joblist.session{job_ind}, job_name);
    fullfilename = [batch_dir '/' filename];
    
    disp(['Creating Batch file for Job ' num2str(job_ind) ': ' fullfilename])
    if exist(fullfilename, 'file')
        if ~strcmpi(overwrite_file, 'a')
            disp(' ')
            disp(['File ' fullfilename ' already exist.'])
            overwrite_file = input('Should it be overwritten? [y, n, a (all)]: ', 's');
            if ~(strcmpi(overwrite_file, 'y') || strcmpi(overwrite_file, 'a'))
                error(['File ' filename ' already exists and should not be overwritten. Solve problem and start again.'])
            end
        end
        delete(fullfilename)
    end
    
    % open single subject file
    fid_single = fopen(fullfilename , 'w');
    
    % ensure that the right shell is used !#/bin/bash
    fprintf(fid_single, '#PBS -S /bin/bash\n');
    % add a comment what this script does
    jobnameline = ['#PBS -N ' job_name '\n'];
    try
        fprintf(fid_single, jobnameline);
    catch ME
        display(ME)
    end
    
    fprintf(fid_single, '#PBS -j oe\n');
    fprintf(fid_single, '#PBS -lnodes=1:ppn=16\n');
    fprintf(fid_single, '#PBS -lwalltime=48:00:00\n');
    
    fprintf(fid_single, '#PBS -o $HOME/PRF/Logs/\n');
    fprintf(fid_single, '#\n');
    
    fprintf(fid_single,['cp -r $HOME/PRF/Data/ses-' joblist.session{job_ind} '* "$TMPDIR"\n']);
    fprintf(fid_single, 'cp -r $HOME/PRF/Code/BashScripts "$TMPDIR"\n');
    fprintf(fid_single, 'cp -r $HOME/PRF/Code/analyzePRF "$TMPDIR"\n');
    fprintf(fid_single, 'cp -r $HOME/PRF/Code/NIfTI "$TMPDIR"\n');
    % get the command to start the job
    % this command will be saved in the job script
    
    fprintf(fid_single,'cd "$TMPDIR"\n');
    fprintf(fid_single,['chmod +x ' execute_matlab_process_sh '\n']);
    line = sprintf('%s %s %s %s %s', execute_matlab_process_sh, parallel_fun, ...
        joblist.session{job_ind}, log_file_dir, parallel_fun_dir);
    fprintf(fid_single, '%s\n', line);
    
    % finally: pass exit status of execute_matlab_process.sh to LISA
    fprintf(fid_single, 'exit $?\n');
    fclose(fid_single);
    
    % [CK] Don't think I need this... ======
    % maxcounter = 0; nrbashscr = 1;
    
    disp(['Adding ' fullfilename ' to original batch file.'])
    % add this script to the list off all scripts
    % be careful to use the linux path
    fprintf(fid_commit_ori, '# Job %s\n', joblist.session{job_ind});
    
    % [CK] Don't think I need this... ======
    %waittime = job_ind*job_ind2*job_ind3*job_ind4; %Waiting time in seconds
    %fprintf(fid_commit_ori,'echo "Waiting for %d seconds to avoid race condition..."\n',waittime);
    %fprintf(fid_commit_ori,'sleep %ds\n',waittime);
    %fprintf(fid_commit_ori,'echo "Waiting done, starting process..."');
    
    line = sprintf('%s %s', 'qsub ', fullfilename);
    fprintf(fid_commit_ori, '%s\n\n', line);
    
    % [CK] Don't think I need this... ======
    %         %% add the produced single job batch file to the main batch file
    %         disp(['Adding ' fullfilename ' to boot batch file.'])
    %         % add this script to the list off all scripts
    %         % be careful to use the linux path
    %         fprintf(fid_commit_all, '# Job %s\n', joblist.session{job_ind});
    %         %waittime = job_ind*job_ind2*job_ind3*job_ind4; %Waiting time in seconds
    %         %fprintf(fid_commit_ori,'echo "Waiting for %d seconds to avoid race condition..."\n',waittime);
    %         %fprintf(fid_commit_ori,'sleep %ds\n',waittime);
    %         %fprintf(fid_commit_ori,'echo "Waiting done, starting process..."');
    %
    %         line = sprintf('%s %s', 'qsub ', fullfilename);
    %         fprintf(fid_commit_all, '%s\n\n', line);
    % =======================================
end
fclose(fid_commit_all);
fclose(fid_commit_ori);
end
