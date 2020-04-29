tool_basepath = '~/Dropbox/MATLAB_NONGIT/TOOLBOX';
addpath(genpath(fullfile(tool_basepath, 'NIfTI')));

homefld = pwd;

MONKEYS={'danny','eddy'};

for m=1:length(MONKEYS)
    MONKEY=MONKEYS{m};

    fprintf(['Doing things for ' MONKEY '\n']);
    ResFld =  ['~/Documents/MRI_ANALYSIS/pRF-NHP/Results/' MONKEY '/'];
    % ResFld =  ['~/Documents/MRI_ANALYSIS/pRF-NHP/Results/LISA/' MONKEY '/'];
    Ref = ['~/Documents/MRI_ANALYSIS/pRF-NHP/Results/Reference/' MONKEY '/'....
        'output_files/Func/nii/func_brain.nii'];

    cd(ResFld);
    fld = dir('2*');

    % load reference nifti
    nii_ref = load_nii(Ref); ref_hdr = nii_ref.hdr;
        
    for s=1:length(fld)
        cd(fullfile(fld(s).folder, fld(s).name));
        matfile = ls('*mat');
        fprintf(['Loading ' matfile '\n']);
        R = load(matfile(1:end-4));
        R.result.ref_hdr = ref_hdr;
        % save the result ----
        result=R.result;
        save(matfile(1:end-1),'result','-v7.3');
        % also save as nifti files
        % angle ---
        fprintf('Angles ');
        nii = make_nii(result.ang,[1 1 1],[],16,...
            'pRF fit: Angles (deg)');
        nii.hdr.hist = ref_hdr.hist;
        save_nii(nii, fullfile(fld(s).folder, fld(s).name, ...
            [matfile(5:17) '_ang.nii']));
        gzip(fullfile(fld(s).folder, fld(s).name, ...
            [matfile(5:17) '_ang.nii']));
        delete(fullfile(fld(s).folder, fld(s).name, ...
            [matfile(5:17) '_ang.nii']));
        % ecc ---
        fprintf('Ecc ');
        nii = make_nii(result.ecc,[1 1 1],[],16,...
            'pRF fit: Eccentricity (pix)');
        nii.hdr.hist = ref_hdr.hist;
        save_nii(nii, fullfile(fld(s).folder, fld(s).name, ...
            [matfile(5:17) '_ecc.nii']));
        gzip(fullfile(fld(s).folder, fld(s).name, ...
            [matfile(5:17) '_ecc.nii']));
        delete(fullfile(fld(s).folder, fld(s).name, ...
            [matfile(5:17) '_ecc.nii']));
        % size ---
        fprintf('Size ');
        nii = make_nii(result.rfsize,[1 1 1],[],16,...
            'pRF fit: RF size (pix)');
        nii.hdr.hist = ref_hdr.hist;
        save_nii(nii, fullfile(fld(s).folder, fld(s).name, ...
            [matfile(5:17) '_rfsize.nii']));
        gzip(fullfile(fld(s).folder, fld(s).name, ...
            [matfile(5:17) '_rfsize.nii']));
        delete(fullfile(fld(s).folder, fld(s).name, ...
            [matfile(5:17) '_rfsize.nii']));
        % R^2 Goodness of fit ---
        fprintf('R2 ');
        nii = make_nii(result.R2,[1 1 1],[],16,...
            'pRF fit: R2 Goodnes off fit');
        nii.hdr.hist = ref_hdr.hist;
        save_nii(nii, fullfile(fld(s).folder, fld(s).name, ...
            [matfile(5:17) '_R2.nii']));
        gzip(fullfile(fld(s).folder, fld(s).name, ...
            [matfile(5:17) '_R2.nii']));
        delete(fullfile(fld(s).folder, fld(s).name, ...
            [matfile(5:17) '_R2.nii']));
        fprintf('>>> DONE\n');
    end
end

cd(homefld);