#!/bin/bash
#SBATCH -N 1 --ntasks-per-node=16
#SBATCH -t 48:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=p.c.klink@gmail.com

source ~/.bash_profile
source ~/.bashrc
umask u+rwx,g+rwx

echo job id $SLURM_JOBID
echo job name $SLURM_JOB_NAME
echo submitted by $SLURM_JOB_ACCOUNT
echo from $SLURM_SUBMIT_DIR
echo the allocated nodes are: $SLURM_JOB_NODELIST

# INFO: FitPRF_avg_dhrf_AllSessions-avg-odd_29:42

mkdir -p $TMPDIR/PRF
mkdir -p $TMPDIR/PRF/Logs/
cp -r $HOME/PRF/Data/avg/danny/AllSessions-avg-odd* $TMPDIR/PRF
cp -r $HOME/PRF/Data/mask/danny/* $TMPDIR/PRF
cp -r $HOME/PRF/Data/refhdr/danny* $TMPDIR/PRF
cp -r $HOME/PRF/Code/* $TMPDIR/PRF
cd $TMPDIR/PRF

chmod +x $TMPDIR/PRF/BashScripts/pRF_run_analyzePRF_LISA_avg_dhrf.sh

$TMPDIR/PRF/BashScripts/pRF_run_analyzePRF_LISA_avg_dhrf.sh \
	pRF_FitModel_LISA_avg danny AllSessions-avg-odd 29:42 defaultHRF [] \
	/home/pcklink/PRF/Logs/ \
	$TMPDIR/PRF/ |& tee $TMPDIR/PRF/Logs/Log_danny_AllSessions-avg-odd_29:42_defaultHRF.txt

cp $TMPDIR/PRF/Logs/Log_danny_AllSessions-avg-odd_29:42_defaultHRF.txt $HOME/PRF/Logs/
exit $?
