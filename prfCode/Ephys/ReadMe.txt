====================================
Workflow EPHYS pRF's
====================================

These pre-fitting scripts are located in the <..>/prfCode/Ephys folder

1	Xing gave me the preprocessed data
	- See EphysNotes.txt for some specs

2	ck_Run.m
	- Central script to coordinate data handling
	- Uses:
		
	2a	ck_Load.m
		- Loads data from <DATA@SERVER>/PRF_EPHYS/
			- Data_raw / Data_preproc / Log_pRF / Channelmaps
		- Saves in Data_proc

3	Copy <...>/Data/ephys/ to LISA <...>/PRF/Data/ephys
	Run fits on LISA
	Copy fit results from LISA <...>/PRF/Results/ back to local FitResults/ephys/<subject>/<model>

4	Add the classic RF mapping data
	- Receive mat files from Xing
	- Copy them to FitResults/ephys/<subject>/classicRF

5	Continue analysis in <..>/prfCode/PostFit