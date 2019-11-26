module load matlab
echo "mcc -m $HOME/PRF/Code/pRF_FitModel_LISA_cv.m -a $HOME/PRF/Code/analyzePRF -a $HOME/PRF/Code/HRF -a $HOME/PRF/Code/NIfTI" | matlab -nodisplay
#echo "mcc -m $HOME/PRF/Code/pRF_FitModel_LISA_ephys.m -a $HOME/PRF/Code/analyzePRF" | matlab -nodisplay
#echo "mcc -m $HOME/PRF/Code/pRF_FitModel_LISA_ephys_neggain.m -a $HOME/PRF/Code/analyzePRF" | matlab -nodisplay
module unload matlab