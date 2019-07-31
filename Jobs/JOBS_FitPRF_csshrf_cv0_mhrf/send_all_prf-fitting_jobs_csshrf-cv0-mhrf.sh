#!/bin/bash
#
# This bash-script submits all jobs to the server, instead of running them locally.
# If you want to submit only some jobs to the server,simply add a "#" in front of 
#the ones you like to ommit and execute the script then.
#

mkdir -p $HOME/PRF/Logs/slurm
cd $HOME/PRF/Logs/slurm
chmod +x $HOME/PRF/Code/Jobs/*

sbatch  $HOME/PRF/Code/Jobs/run_job_Ses-AllSessions-avg-cv0_1_FitPRF_csshrf_cv0_mhrf_danny.sh
sbatch  $HOME/PRF/Code/Jobs/run_job_Ses-AllSessions-avg-cv0_2_FitPRF_csshrf_cv0_mhrf_danny.sh
sbatch  $HOME/PRF/Code/Jobs/run_job_Ses-AllSessions-avg-cv0_3_FitPRF_csshrf_cv0_mhrf_danny.sh
sbatch  $HOME/PRF/Code/Jobs/run_job_Ses-AllSessions-avg-cv0_4_FitPRF_csshrf_cv0_mhrf_danny.sh

sbatch  $HOME/PRF/Code/Jobs/run_job_Ses-AllSessions-avg-cv0_1_FitPRF_csshrf_cv0_mhrf_eddy.sh
sbatch  $HOME/PRF/Code/Jobs/run_job_Ses-AllSessions-avg-cv0_2_FitPRF_csshrf_cv0_mhrf_eddy.sh
sbatch  $HOME/PRF/Code/Jobs/run_job_Ses-AllSessions-avg-cv0_3_FitPRF_csshrf_cv0_mhrf_eddy.sh
sbatch  $HOME/PRF/Code/Jobs/run_job_Ses-AllSessions-avg-cv0_4_FitPRF_csshrf_cv0_mhrf_eddy.sh


