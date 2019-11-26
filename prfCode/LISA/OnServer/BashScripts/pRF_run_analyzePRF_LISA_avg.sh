#! /bin/bash

# compiles & calls the compiled matlab file specified as first parameter, 
# passes the 2nd to 6th parameters as input

# PARAMETERS ====================
echo 'Matlab file name: '$1
MATLAB_SCRIPT=$1

echo 'Monkey name: '$2
MONKEY=$2

echo 'Session: '$3
SESS=$3

echo 'Slices: '$4
SLICES=$4
# create folder label
SLICELABEL=${SLICES/:/-}

echo 'HRF: '$5
HRF=$5

echo 'Number of workers: '$6
NWORKERS=$6

echo 'Modeltype: '$7
MODTYPE=$7

echo 'Crossvalidation mode: '$8
XVAL=$8

echo 'Results folder: '$9
RESFLD=$9

echo 'Path to Logfile directory: '$10
LOGFILE_DIR=$10

echo 'Path to matlab file: '$11
MATLAB_PATH=$11


# check the number of arguments we need
NARG=11 
if [ $# -ne $NARG ];then
    echo "ERROR: You must exactly provide $NARG arguments!" 1>&2
    exit 1
fi

LOC=`pwd`
cd $HOME 
# ensure that you are in your home directory, to avoid matlab path problems 
# (by default, the path information is saved in the home directory, so starting 
# matlab from there uses the normal path of the user)

HOSTNAME=`hostname`
MYSELF=`whoami`
DATE=`date`

LOGFILE=$LOGFILE_DIR/log_job_$MONKEY_$SESS_$SLICES

echo
echo "Running job as $MYSELF on $HOSTNAME, $DATE" 2>&1
echo

# creating a new virtual screen using xvfb-run and calling then matlab, 
# including a 'mini-program' to ensure all paths are added correctly
# " , \ " at end of line is necessary, because this whole little program 
# must be transmitted as one line to matlab.
echo

# go to the scratch PRF folder
cd $TMPDIR/PRF
# copy the matlab script that needs to be compiled from Code/BashScripts to scratch/PRF
cp $HOME/PRF/Code/$MATLAB_SCRIPT* $TMPDIR/PRF 

# compile the matlab script (only need to do this once in theory)
#module load matlab
#echo "mcc -m $TMPDIR/PRF/$MATLAB_SCRIPT.m -a $TMPDIR/PRF/analyzePRF -a $TMPDIR/PRF/HRF -a $TMPDIR/PRF/NIfTI" | matlab -nodisplay
#module unload matlab
#wait

# make it executable
chmod +x $MATLAB_SCRIPT*

# load the mcr module (runtime environment)
module load mcr
# execute the compiled function with the specified inputs
./$MATLAB_SCRIPT $MONKEY $SESS $SLICES $HRF $NWORKERS $MODTYPE $XVAL $RESFLD
# unload the runtime environment
module unload mcr
# change back to the starting directory
cd $LOC

echo
echo "Job finished with Exit code $EXIT, returning"
echo "Copy data from scratch back to home " 2>&1
echo

mkdir -p $HOME/PRF/Results/$MONKEY/$RESFLD
cp -r $TMPDIR/PRF/Results/$MONKEY/$RESFLD/Slices* $HOME/PRF/Results/$MONKEY/$RESFLD/

exit $EXIT
