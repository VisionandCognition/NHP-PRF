#! /bin/bash

# compiles & calls the compiled matlab file specified as first parameter, 
# passes the 2nd and 3rd parameters as input

# PARAMETERS ====================
echo 'Matlab file name: '$1
MATLAB_SCRIPT=$1

echo 'Monkey name: '$2
MONKEY=$2

echo 'Session : '$3
SESS=$3

echo 'Path to Logfile directory: '$4
LOGFILE_DIR=$4

echo 'Path to matlab file: '$5
MATLAB_PATH=$5


# check the number of arguments we need
NARG=5
if [ $# -ne $NARG ];then
    echo "ERROR: You must exactly provide $NARG arguments!" 1>&2
    exit 1
fi

LOC=`pwd`
cd "$HOME" 
# ensure that you are in your home directory, to avoid matlab path problems 
# (by default, the path information is saved in the home directory, so starting 
# matlab from there uses the normal path of the user)

HOSTNAME=`hostname`
MYSELF=`whoami`
DATE=`date`

LOGFILE=$LOGFILE_DIR/log_job_$MONKEY$SESS

echo
echo "Running job as $MYSELF on $HOSTNAME, $DATE" 2>&1
echo


# creating a new virtual screen using xvfb-run and calling then matlab, 
# including a 'mini-program' to ensure all paths are added correctly
# " , \ " at end of line is necessary, because this whole little program 
# must be transmitted as one line to matlab.
echo

# go to the scratch PRF folder
cd "$TMPDIR"/PRF
# copy the matlab script that needs to be compiled from Code/BashScripts to scratch/PRF
cp $HOME/PRF/Code/$MATLAB_SCRIPT "$TMPDIR"/PRF 

# compile the matlab script
module load matlab
mcc -m "$TMPDIR"/PRF/$MATLAB_SCRIPT.m -a "$TMPDIR"/PRF/analyzePRF -a "$TMPDIR"/PRF/NIfTI
module unload matlab

# make it executable
chmod +x $MATLAB_SCRIPT

# load the mcr module (runtime environment)
module load mcr
# execute the compiled function with the specified inputs
./$MATLAB_SCRIPT $MONKEY $SESS
# unload the runtime environment
module unload mcr
# change back to the starting directory
cd $LOC

echo
echo "Job finished with Exit code $EXIT, returning"

  echo "Copy data from scratch back to home " 2>&1
echo

mkdir $HOME/PRF/Results/$MONKEY
cp -r "$TMPDIR"/PRF/Results/%MONKEY/$SESS $HOME/PRF/Results/$MONKEY


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
