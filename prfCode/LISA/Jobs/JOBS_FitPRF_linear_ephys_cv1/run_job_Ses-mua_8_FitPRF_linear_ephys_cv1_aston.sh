#!/bin/bash
#SBATCH -N 1 --ntasks-per-node=16
#SBATCH -t 05:00:00
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

# INFO: FitPRF_linear_ephys_cv1_mua_8

mkdir -p $TMPDIR/PRF
mkdir -p $TMPDIR/PRF/Logs/
cp -r $HOME/PRF/Code/* $TMPDIR/PRF
cd $TMPDIR/PRF

chmod +x $TMPDIR/PRF/BashScripts/pRF_run_analyzePRF_LISA_ephys.sh

$TMPDIR/PRF/BashScripts/pRF_run_analyzePRF_LISA_ephys.sh \
	pRF_FitModel_LISA_multephys aston mua 8 [] \
	linear_ephys 1 linear_ephys_cv1 /home/pcklink/PRF/Logs/ $TMPDIR/PRF/ \
	 |& tee $TMPDIR/PRF/Logs/Log_aston_mua_8_linear_ephys_xval1.txt

cp $TMPDIR/PRF/Logs/Log_aston_mua_8_linear_ephys_xval1.txt $HOME/PRF/Logs/
exit $?
