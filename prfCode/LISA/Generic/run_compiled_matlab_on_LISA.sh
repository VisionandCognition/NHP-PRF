#! /bin/bash

# calls the compiled matlab file specified as first parameter, passing the second parameter, e.g.
# example_script(5)
#
# PARAMETERS
echo $1': matlab file name (e.g. "example_script", NOT example_script.m)'
MATLAB_SCRIPT=$1

echo $2': Job_nr'
JOB_NR1=$2
cd ..
echo $3': Job_nr2'
JOB_NR2=$3

echo $4': Job_nr3'
JOB_NR3=$4

echo $5': Job_nr4'
JOB_NR4=$5

echo $6': Path to Logfile directory'
LOGFILE_DIR=$6

echo $7': Path to matlab file (e.g. /analysis/tom/example_dir/)'
MATLAB_PATH=$7


# check the number of arguments we need
NARG=7
if [ $# -ne $NARG ];then
    echo "ERROR: You must exactly provide $NARG arguments (see execute_matlab_process.sh for details)!" 1>&2
    exit 1
fi

# 2nd argument must be an integer
echo $(( $JOB_NR1 + 1 )) >& /dev/null
if [ $? -ne 0 ];then
    echo "ERROR: Second argument must be a job number (integer)!" 1>&2
    exit 1
fi

%

LOC=`pwd`
cd "$HOME" # ensure that you are in your home directory, to avoid matlab path problems (by default, the path information is saved in the home directory, so staring matlab from there uses the normal path of the user)

HOSTNAME=`hostname`
MYSELF=`whoami`
DATE=`date`

LOGFILE=$LOGFILE_DIR/log_job$JOB_NR1$JOB_NR2$JOB_NR3$JOB_NR4

	echo
echo "Running job as $MYSELF on $HOSTNAME, $DATE" 2>&1
echo



# creating a new virtual screen using xvfb-run and calling then matlab, including a 'mini-program' to ensure all paths are added correctly
# " , \ " at end of line is necessary, because this whole little program must be transmitted as one line to matlab.
echo

cd "$TMPDIR"
cp $HOME/MVPA_Scripts/$MATLAB_SCRIPT "$TMPDIR" # Copy all scripts 

mkdir TMPResults
mkdir TMPResults/TMPMatlab

cp $HOME/TMPResults/TMPMatlab/$JOB_NR1"_"$JOB_NR3"_tmpfile.mat" "$TMPDIR"/TMPResults/TMPMatlab/ 


#module load matlab

#mcc -m "$TMPDIR"/MVPA_Scripts/$MATLAB_SCRIPT.m -a "$TMPDIR"/MVPA_Scripts/*.m
#module unload matlab

chmod +x $MATLAB_SCRIPT


module load mcr
./$MATLAB_SCRIPT $JOB_NR1 $JOB_NR2 $JOB_NR3 $JOB_NR4

module unload mcr
# change back to the directory we came from
cd $LOC

echo
echo "Job finished with Exit code $EXIT, returning"


  echo "Copy data from scratch back to home " 2>&1
echo

#cp -r "$TMPDIR"/TMPData/*  $HOME/TMPData/
cp -r "$TMPDIR"/TMPResults/* $HOME/TMPResults/
if [ $EXIT -ne 0 ] ; then
  echo "WARNING: Job finished with error code $EXIT. 
Check the logfile why.

If you have no clue why, and this happend the first time , just try to 
restart this job (copy the lines from to the terminal 
../batch_dir/send_condor_all_jobs.sh). Dont call send_condor_all_jobs.sh 
again, because this will restart all jobs.

If it crashs again, check if all machines are busy, and if so, wait a bit, 
and restart it again.

If no normal logfile was created, the path to matlab that you provided 
might be wrong ($MATLAB_EXE). If this happens often, the target 
machine might be out of memory.

If your jobs always crash on the same machine, copy & rename
condor_generate_submit.sh and create your own version in which you exclude 
especially this machine (see condor_generate_submit_kai.sh as an example),
then modify create_parallel.m to use this script.

Job ran as $MYSELF on $HOSTNAME, $DATE" >> $LOGFILE.err
fi

# and finally return the exit status we got from xvfb
exit $EXIT
