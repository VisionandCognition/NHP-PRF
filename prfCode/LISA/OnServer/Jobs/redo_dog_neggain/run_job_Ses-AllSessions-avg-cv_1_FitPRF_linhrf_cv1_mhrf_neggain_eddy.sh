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
module load pre2019

# INFO: FitPRF_linhrf_cv1_mhrf_neggain_AllSessions-avg-cv_01:15

mkdir -p $TMPDIR/PRF
mkdir -p $TMPDIR/PRF/Logs/
cp -r $HOME/PRF/Data/cv/eddy/AllSessions-avg-cv* $TMPDIR/PRF
cp -r $HOME/PRF/Data/mask/eddy/* $TMPDIR/PRF
cp -r $HOME/PRF/Data/refhdr/eddy* $TMPDIR/PRF
cp -r $HOME/PRF/Code/* $TMPDIR/PRF
cd $TMPDIR/PRF

chmod +x $TMPDIR/PRF/BashScripts/pRF_run_analyzePRF_LISA_avg.sh

$TMPDIR/PRF/BashScripts/pRF_run_analyzePRF_LISA_avg.sh \
	pRF_FitModel_LISA_cv_neggain eddy AllSessions-avg-cv 01:15 HRF_monkey [] \
	linear_hrf 1 linhrf_cv1_mhrf_neggain /home/pcklink/PRF/Logs/ $TMPDIR/PRF/ \
	 |& tee $TMPDIR/PRF/Logs/Log_eddy_AllSessions-avg-cv_01:15_HRF_monkey_linear_hrf_xval1.txt

cp $TMPDIR/PRF/Logs/Log_eddy_AllSessions-avg-cv_01:15_HRF_monkey_linear_hrf_xval1.txt $HOME/PRF/Logs/
exit $?
