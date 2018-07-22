#PBS -S /bin/bash
#PBS -N FitPRF_PerSession_20160804_31:45
#PBS -j oe
#PBS -lnodes=1:ppn=16:mem64gb
#PBS -lwalltime=48:00:00
#PBS -o $HOME/PRF/Logs/
#

echo "Job $PBS_JOBID started at `date`. Subject: eddy, Session: 20160804, Slices: 31:45" | mail $USER -s "Job $PBS_JOBID"

mkdir $TMPDIR/PRF
cp -r $HOME/PRF/Data/us_reg/eddy/ses-20160804* $TMPDIR/PRF
cp -r $HOME/PRF/Data/mask/eddy/* $TMPDIR/PRF
cp -r $HOME/PRF/Code/* $TMPDIR/PRF
cd $TMPDIR/PRF

chmod +x $TMPDIR/PRF/BashScripts/pRF_run_CompiledMatlab_LISA.sh

$TMPDIR/PRF/BashScripts/pRF_run_CompiledMatlab_LISA.sh \
	pRF_FitModel_LISA eddy 20160804 31:45 [] \
	/home/pcklink/PRF/Logs/ \
	$TMPDIR/PRF/


echo "Job $PBS_JOBID ended at `date`. Subject: eddy, Session: 20160804, Slices: 31:45" | mail $USER -s "Job $PBS_JOBID"

exit $?
