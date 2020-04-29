====================================
Workflow POSTFIT pRF's (matlab)
====================================

These post-fitting scripts are located in the <..>/prfCode/PostFit/matlab folder

1	BATCH_ck_Combine_MRI_SliceChunks.m
	- Uses ck_Combine_MRI_SliceChunks.m to put chunks back together
	- Takes slice chunks from FitResults/MRI/<subject>/<model>
	- Combine them
	- Save volumes in same folder

2	ck_CollectAllFittingResults.m
	- ck_GetMRI_pRF.m
		- combines results from 2 animals
		- Gets data from <..>/FitResults/MRI/<subject> saves data in <..>/FitResults/MRI/Combined
		- Adds D99 atlasses
			- ck_GetAtlasTable.m gets the atlas info
	- ck_GetEphys_pRF.m
		- combines results from 2 animals
		- Gets and saves data in <..>/FitResults/ephys
 	- Combines MRI and Ephys
 	- Saves as AllFits.cv1.mat in <..>/FitResults/MultiModal

3	Run ck_BasicProc_AllFits_cv1.m
	- Takes data from <..>/FitResults/MultiModal	
	- Re-structures
	- Saves back as tables and structures in <..>/FitResults/MultiModal
		
4	ck_StatsAndPlots.m
	- Long script that creates plots and runs statistics etc
	- Uses: 
		- ck_GetROIidx.m
			- get an ROI index back when querying with label
		- ck_xmod_stats.m
			- Runs the LMM for crossmodal statistics
