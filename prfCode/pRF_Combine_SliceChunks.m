function pRF_Combine_SliceChunks(Monkey,Session)
% For prf analysis that is run on LISA in separate slice chunks, this
% function recombines the results in volumes

tool_basepath = '/Users/chris/Documents/MATLAB/TOOLBOX';
addpath(genpath(fullfile(tool_basepath, 'NIfTI')));

homefld = pwd;
ResFld =  ['~/Documents/MRI_ANALYSIS/NHP-pRF/Results/LISA/' Monkey '/' Session];
cd(ResFld);
fld = dir('Slices_*');
for s=1:length(fld)
    cd(fullfile(fld(s).folder, fld(s).name));
    matfile = ls('*mat');
    fprintf(['Loading ' matfile 'Chunk ' num2str(s) '\n']);
    R = load(matfile(1:end-4));
    
    if s==1
        result = R.result;
    else % concat
        result.ang      = cat(3,result.ang,R.result.ang);
        result.ecc      = cat(3,result.ecc,R.result.ecc);
        result.expt     = cat(3,result.expt,R.result.expt);
        result.rfsize   = cat(3,result.rfsize,R.result.rfsize);
        result.R2       = cat(3,result.R2,R.result.R2);
        result.gain     = cat(3,result.gain,R.result.gain);
        result.resnorms = cat(3,result.resnorms,R.result.resnorms);
        result.numiters = cat(3,result.numiters,R.result.numiters);
        result.meanvol  = cat(3,result.meanvol,R.result.meanvol);
    end
    clear R
end
cd(ResFld);

% save the result ----
fprintf('Merging and saving the result: mat-file, ');
save(fullfile(ResFld,['pRF_Sess-' Session]),'result','-v7.3');

% also save as nifti files
% angle ---
fprintf('Angles, ');
nii = make_nii(result.ang,[1 1 1],[],16,...
    'pRF fit: Angles (deg)');
nii.hdr.hist = result.hdr_ref.hist;
save_nii(nii, fullfile(ResFld, ['Sess-' Session '_ang.nii']));
gzip(fullfile(ResFld, ['Sess-' Session '_ang.nii']));
delete(fullfile(ResFld, ['Sess-' Session '_ang.nii']));
% ecc ---
fprintf('Ecc, ');
nii = make_nii(result.ecc,[1 1 1],[],16,...
    'pRF fit: Eccentricity (pix)');
nii.hdr.hist = result.hdr_ref.hist;
save_nii(nii, fullfile(ResFld, ['Sess-' Session '_ecc.nii']));
gzip(fullfile(ResFld, ['Sess-' Session '_ecc.nii']));
delete(fullfile(ResFld, ['Sess-' Session '_ecc.nii']));
% size ---
fprintf('Size, ');
nii = make_nii(result.rfsize,[1 1 1],[],16,...
    'pRF fit: RF size (pix)');
nii.hdr.hist = result.hdr_ref.hist;
save_nii(nii, fullfile(ResFld, ['Sess-' Session '_rfsize.nii']));
gzip(fullfile(ResFld, ['Sess-' Session '_rfsize.nii']));
delete(fullfile(ResFld, ['Sess-' Session '_rfsize.nii']));
% R^2 Goodness of fit ---
fprintf('R2 ');
nii = make_nii(result.R2,[1 1 1],[],16,...
    'pRF fit: R2 Goodnes off fit');
nii.hdr.hist = result.hdr_ref.hist;
save_nii(nii, fullfile(ResFld, ['Sess-' Session '_R2.nii']));
gzip(fullfile(ResFld, ['Sess-' Session '_R2.nii']));
delete(fullfile(ResFld, ['Sess-' Session '_R2.nii']));

fprintf('>> Done!\n');
cd(homefld)
