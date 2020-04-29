====================================
Workflow MRI pRF's
====================================

These pre-fitting scripts are located in the <..>/prfCode/MRI folder

1	pRF_prepdata_avg.m (batch: BATCH_pRF_PadAndUpsample.m)
	- Takes pRF_PrepDatalist_MONKEY.m as info on what preprocessed fMRI data to load
	- Creates timeseries of equal 230 volumes length by inserting NaNs
	- Saves as ses-XXXX-230vols.mat 
	- In <...>/Data/MRI/us-padded/<subject>/

2	pRF_sessavg_BOLD_tseries.m (batch: BATCH_pRF_avg_BOLD_tseries.m)
	- Takes the ses-XXXX-230vols.mat from <...>/Data/MRI/us-padded/<subject>/
	- Creates an average timecourse per session
	- Saves as medianBOLD_sess-xxxxxx.mat
	- In <...>/Data/MRI/us-padded/<subject>/

3	pRF_AvgSessions.m / pRF_AvgSessions_odd.m / pRF_AvgSessions_even.m
	- Takes the medianBOLD_sess-xxxxxx.mat from <...>/Data/MRI/us-padded/<subject>/
	- Averages across sessions
	- Saves as AllSessions-avg.mat / AllSessions-avg-odd.mat / AllSessions-avg-even.mat
	- Copies these files to <...>/Data/MRI/avg/<subject>/

4	pRFprepCrossVal.m
	- Takes AllSessions-avg-odd.mat & AllSessions-avg-even.mat from <...>/Data/MRI/avg/<subject>/
	- Restructures
	- Saves as AllSessions-avg-cv.mat
	- In <...>/Data/MRI/cv/<subject>/

5	Copy <...>/Data/MRI/cv to LISA <...>/PRF/Data/MRI/cv
	Run fits on LISA
	Copy fit results from LISA <...>/PRF/Results/ back to local FitResults/MRI/<subject>/<model>

6	Continue analysis in <..>/prfCode/PostFit