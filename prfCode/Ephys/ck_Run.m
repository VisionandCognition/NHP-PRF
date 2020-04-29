%% ck_Run
% use this file to run subfunctions
clc;
addpath(genpath('../../prfCode'));
addpath(fullfile('/home','chris','Documents','MATLAB','GENERAL_FUNCTIONS'));

%% SETTINGS ###############################################################
Subj = {'Aston','Lick'}; % Lick / Aston
Sess = {'20181004_B1','20180807_B2'}; % 20180807_B2 / 20181004_B1

for s=2:length(Subj)
    subj = Subj{s};
    sess = Sess{s};
    
    % MORE SETTINGS #######################################################
    Do.Load.Any         = true; 
    Do.Load.DigChan     = true; % new files
    Do.Load.MUA         = true; % new files
    Do.Load.LFP         = true; % new files
    Do.Load.Behavior    = true; % new files

    Do.Load.ProcMUA     = false;
    Do.Load.ProcLFP     = false;
    Do.Load.ProcBEH     = false;

    Do.SyncTimes        = true;
    Do.SaveUncut        = true;

    Do.SaveMUA_perArray = true;
    Do.SaveLFP_perArray = true;

    % Switched to AnalyzePRF, do not use this -----
    Do.CreatePrediction = false; 

    Do.FitPRF           = false;
    Do.FitMUA           = false; % MUA data
    Do.FitLFP           = false; % LFP data (split by freq band)

    Do.PlotPRF_MUA      = false;
    Do.PlotPRF_LFP      = false;
    % ---------------------------------------------
    

    %% load data ==========================================================
    %This preprocess the data so it can be used in fits.
    if Do.Load.Any
        ck_Load(subj, sess, Do)
    end



    %% ####################################################################
    % This is the old way of doing it. Creates a large range of predicted
    % pRFs. The fitting will then find the best pRF per channel/signal-type
    % NB1! Switched to using the adapted AnalyzePRF toolbox.
    % >> This is no longer used <<

    % create pRF prediction ==============================================
    % NB2! Should be done in parallel and will take days.
    if Do.CreatePrediction
        ck_CreatePRFPrediction(subj, sess) 
    end

    % do the fitting procedure ===========================================

    if Do.FitPRF
        ck_FitPRF(subj, sess, Do)
    end

    % plot some results ==================================================
    if Do.PlotPRF_MUA
        ck_PlotPRF_MUA(subj,sess)
    end
    if Do.PlotPRF_LFP
        ck_PlotPRF_LFP(subj,sess)
    end

    % #####################################################################

end
rmpath(genpath('../../prfCode'));
