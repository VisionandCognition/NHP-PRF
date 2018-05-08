#PBS -S /bin/bash
#PBS -N FitPRF_PerSession
#PBS -j oe
#PBS -lnodes=1:ppn=16
#PBS -lwalltime=48:00:00
#PBS -o $HOME/PRF/Logs/
#
cp -r $HOME/PRF/Data/us_reg/ses-20180117* "$TMPDIR"/PRF
cp -r $HOME/PRF/Code/BashScripts "$TMPDIR"/PRF
cp -r $HOME/PRF/Code/analyzePRF "$TMPDIR"/PRF
cp -r $HOME/PRF/Code/NIfTI "$TMPDIR"/PRF

cd "$TMPDIR/PRF"

chmod +x /media/DOCUMENTS/DOCUMENTS/MRI_ANALYSIS/analyzePRF/scratch/PRF/BashScripts/run_compiled_matlab_on_LISA.sh

/media/DOCUMENTS/DOCUMENTS/MRI_ANALYSIS/analyzePRF/scratch/PRF/BashScripts/run_compiled_matlab_on_LISA.sh pRF_FitModel_LISA ses-20180117 /media/DOCUMENTS/DOCUMENTS/MRI_ANALYSIS/analyzePRF/PRF/FitPRF_PerSession/logs/ /media/DOCUMENTS/DOCUMENTS/MRI_ANALYSIS/analyzePRF/scratch/PRF/

exit $?
