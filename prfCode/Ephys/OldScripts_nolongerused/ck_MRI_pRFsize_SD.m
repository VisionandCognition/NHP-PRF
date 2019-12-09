tool_basepath = '~/Dropbox/MATLAB_NONGIT/TOOLBOX';
addpath(genpath(fullfile(tool_basepath, 'NIfTI')));

s='/Users/chris/Dropbox/CURRENT_PROJECTS/NHP_MRI/Projects/pRF/FitResults/MRI/org';

m='danny';
f={'20171116','20171129_1','20171129_2','20171207','20171214','20171214_1',...
    '20171214_2','20171220','20180117','20180117_1','20180117_2','20180124',...
    '20180125_1','20180125_2','20180131','20180201'};

for i=1:length(f)
    cd(fullfile(s,m,f{i}));
    c=ls('*.mat');c=c(1:end-1);
    load(c);
    rfs_sig = result.rfsize.*sqrt(result.expt);
    
    nii = make_nii(rfs_sig,[1 1 1],[],16,...
            'pRF fit: RF size (sigma pix)');
    
    nii.hdr.hist = result.ref_hdr.hist;
    
    save_nii(nii, [c(5:17) '_sd.nii']);
end


m='eddy';
f={'20160721','20160728','20160729','20160803','20160804','20170411',...
    '20170512','20170518','20171129'};

for i=1:length(f)
    cd(fullfile(s,m,f{i}));
    c=ls('*.mat');c=c(1:end-1);
    load(c);
    rfs_sig = result.rfsize.*sqrt(result.expt);
    
    nii = make_nii(rfs_sig,[1 1 1],[],16,...
            'pRF fit: RF size (sigma pix)');
    
    nii.hdr.hist = result.ref_hdr.hist;
    
    save_nii(nii, [c(5:17) '_sd.nii']);
end
