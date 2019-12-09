%% ck_Run
% use this file to run subfunctions
addpath(genpath('../../prfCode'));

%% SETTINGS ###############################################################
Subj = {'Lick','Aston'}; % Lick / Aston
Sess = {'20180807_B2','20181004_B1'}; % 20180807_B2 / 20181004_B1

for s=1:length(Subj)
    % MORE SETTINGS #######################################################
    Do.Load.Any         = true; 
    Do.Load.DigChan     = false; % new files
    Do.Load.MUA         = false; % new files
    Do.Load.LFP         = false; % new files
    Do.Load.Behavior    = false; % new files

    Do.Load.ProcMUA     = true;
    Do.Load.ProcLFP     = false;
    Do.Load.ProcBEH     = true;

    Do.SyncTimes        = true;
    Do.SaveUncut        = false;

    Do.SaveMUA_perArray = true;
    Do.SaveLFP_perArray = false;

    Do.CreatePrediction = false; % Switched to AnalyzePRF, do not use this

    Do.FitPRF           = false;
    Do.FitMUA           = false; % MUA data
    Do.FitLFP           = false; % LFP data (split by freq band)

    Do.PlotPRF_MUA      = false;
    Do.PlotPRF_LFP      = false;

    %% load data ==========================================================
    This preprocess the data so it can be used in fits.
    if Do.Load.Any
        ck_Load(subj, sess, Do)
        %ck_Load_Median(subj, sess, Do)
        %ck_Load_Median_firstavgperbar(subj, sess, Do)
    end

    %% ####################################################################
    % This is the old way of doing it. Creates a large range of predicted
    % pRFs. THe fitting will then find the best pRF per channel/signal-type
    % NB1! Switched to using the adapted AnalyzePRF toolbox.
    % >> This is no longer used. 

    %% create pRF prediction ==============================================
    % NB2! Should be done in parallel and will take days.
    if Do.CreatePrediction
        ck_CreatePRFPrediction(subj, sess) 
    end

    %% do the fitting procedure ===========================================

    if Do.FitPRF
        ck_FitPRF(subj, sess, Do)
    end

    %% plot some results ==================================================
    if Do.PlotPRF_MUA
        ck_PlotPRF_MUA(subj,sess)
    end
    if Do.PlotPRF_LFP
        ck_PlotPRF_LFP(subj,sess)
    end

    % #####################################################################

end
rmpath(genpath('../../prfCode'));