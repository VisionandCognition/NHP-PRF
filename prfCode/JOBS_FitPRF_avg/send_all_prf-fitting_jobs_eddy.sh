#!/bin/bash
#
# This bash-script submits all jobs to the server, instead of running them locally.
# If you want to submit only some jobs to the server,simply add a "#" in front of 
#the ones you like to ommit and execute the script then.
#
sbatch  $HOME/PRF/Code/Jobs/run_job_Ses-medianBOLD_sub-eddy_1_FitPRF_avg.sh
sbatch  $HOME/PRF/Code/Jobs/run_job_Ses-medianBOLD_sub-eddy_2_FitPRF_avg.sh
sbatch  $HOME/PRF/Code/Jobs/run_job_Ses-medianBOLD_sub-eddy_3_FitPRF_avg.sh
sbatch  $HOME/PRF/Code/Jobs/run_job_Ses-medianBOLD_sub-eddy_4_FitPRF_avg.sh

