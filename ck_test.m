%% Initialize =============================================================
% Add nifti reading toolbox
addpath(genpath('/Users/chris/Documents/MATLAB/TOOLBOX/NIfTI'))
% Add Kendrick Kay's pRF analysis toolbox
addpath(genpath('/Users/chris/Documents/MRI_ANALYSIS/analyzePRF'))

% Where's the data
datapath = ['/media/DOCUMENTS/Dropbox/CURRENT_PROJECTS/NHP_fMRI/pRF/Data/' ...
    'pRF_example_dataset/Danny20180117'];

% What run to test this?
run_inc = 2;

% Load stimulus file
s = load([datapath '/func/Behavior/stimulus'],'stimulus');

% Load nifti file
for ri = 1:length(run_inc)
    m(ri).nii=load_nii([datapath '/func/mc/run' ...
        num2str(run_inc(ri),'%0.3d') '_mc.feat/filtered_func_data.nii.gz']);
    % NB1 the x-order (DIM 1) is flipped (but anatomically ok)
    m(ri).run = run_inc(ri);
    m(ri).stimulus = s.stimulus{run_inc(ri)};
end

%% process ================================================================
for ri = 1:length(run_inc)
    stimulus{ri} = m(ri).stimulus; % stimulus
    % test the prf mapping on a smaller data set
    slicenr = 30;
    LR = 30;%25:75;
    AP = 15;%3:20;
    data{ri} = NaN(length(LR)*length(AP)*length(slicenr),size(m(ri).nii.img,4));
    voxlist=[]; for x=LR;for y=AP; voxlist = [voxlist; x,y,slicenr]; end; end
    for v=1:size(voxlist,1)
            data{ri}(v,:) = reshape(m(ri).nii.img(...
                voxlist(v,1),voxlist(v,2),voxlist(v,3),:),[1,size(data{ri},2)]);
    end  
end
tr=2.5;

%% run prf analysis
results = analyzePRF(stimulus,data,tr);
