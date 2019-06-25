#!/bin/bash
#SBATCH -N 1 --ntasks-per-node=16
#SBATCH -t 48:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=p.c.klink@gmail.com
#SBATCH -o $HOME/PRF/Logs/


source ~/.bash_profile
source ~/.bashrc
umask u+rwx,g+rwx
echo job id $SLURM_JOBID
echo job name $SLURM_JOB_NAME
echo submitted by $SLURM_JOB_ACCOUNT
echo from $SLURM_SUBMIT_DIR
echo the allocated nodes are: $SLURM_JOB_NODELIST

# INFO: FitPRF_avg_medianBOLD_sub-danny_01:15

mkdir $TMPDIR/PRF
cp -r $HOME/PRF/Data/avg/danny/ses-medianBOLD_sub-danny* $TMPDIR/PRF
cp -r $HOME/PRF/Data/mask/danny/* $TMPDIR/PRF
cp -r $HOME/PRF/Code/* $TMPDIR/PRF
cd $TMPDIR/PRF

chmod +x $TMPDIR/PRF/BashScripts/pRF_run_analyzePRF_LISA_avg.sh

$TMPDIR/PRF/BashScripts/pRF_run_analyzePRF_LISA_avg.sh \
	pRF_FitModel_LISA_avg danny medianBOLD_sub-danny 01:15  [] \
	/home/pcklink/PRF/Logs/ \
	$TMPDIR/PRF/

exit $?
