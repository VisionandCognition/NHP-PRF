#!/bin/bash
#SBATCH -N 1 --ntasks-per-node=16
#SBATCH -t 24:00:00
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

# INFO: FitPRF_csshrf_cv0_dhrf_AllSessions-avg-cv0_43:56

mkdir -p $TMPDIR/PRF
mkdir -p $TMPDIR/PRF/Logs/
cp -r $HOME/PRF/Data/cv0/danny/AllSessions-avg-cv0* $TMPDIR/PRF
cp -r $HOME/PRF/Data/mask/danny/* $TMPDIR/PRF
cp -r $HOME/PRF/Data/refhdr/danny* $TMPDIR/PRF
cp -r $HOME/PRF/Code/* $TMPDIR/PRF
cd $TMPDIR/PRF

chmod +x $TMPDIR/PRF/BashScripts/pRF_run_analyzePRF_LISA_avg.sh

$TMPDIR/PRF/BashScripts/pRF_run_analyzePRF_LISA_avg.sh \
	pRF_FitModel_LISA_cv0 danny AllSessions-avg-cv0 43:56 defaultHRF [] \
	css_hrf 1 csshrf_cv0_dhrf /home/pcklink/PRF/Logs/ $TMPDIR/PRF/ \
	 |& tee $TMPDIR/PRF/Logs/Log_danny_AllSessions-avg-cv0_43:56_defaultHRF_css_hrf_xval1.txt

cp $TMPDIR/PRF/Logs/Log_danny_AllSessions-avg-cv0_43:56_defaultHRF_css_hrf_xval1.txt $HOME/PRF/Logs/
exit $?
