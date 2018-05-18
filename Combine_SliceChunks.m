function Combine_SliceChunks(Monkey,Session)
% For prf analysis that is run on LISA in separate slice chunks, this
% function recombines the results in volumes

tool_basepath = '/Users/chris/Documents/MATLAB/TOOLBOX';
addpath(genpath(fullfile(tool_basepath, 'NIfTI')));

homefld = pwd;
ResFld =  ['~/Documents/MRI_ANALYSIS/pRF-NHP/Results/LISA/' Monkey '/' Session];
cd(ResFld);
fld = dir('Slices_*');
for s=1:length(fld)
    cd(fullfile(fld(s).folder, fld(s).name));
    matfile = ls('*mat');
    fprintf(['Loading ' matfile '\n']);
    R = load(matfile(1:end-4));
    
    if s==1
        result = R.result;
    else % concat
        cat(3,result.ang,R.result.ang);
        cat(3,result.ecc,R.result.ecc);
        cat(3,result.expt,R.result.expt);
        cat(3,result.rfsize,R.result.rfsize);
        cat(3,result.R2,R.result.R2);
        cat(3,result.gain,R.result.gain);
        cat(3,result.resnorms,R.result.resnorms);
        cat(3,result.numiters,R.result.numiters);
        cat(3,result.meanvol,R.result.meanvol);
    end
    clear R
end
cd(ResFld);

% save the result ----
fprintf('Saving the result: ');
save(fullfile(ResFld,['pRF_Sess-' Session]),'result','-v7.3');

% also save as nifti files
% angle ---
fprintf('Angles ');
nii = make_nii(result.ang,[1 1 1],[],[],...
    'pRF fit: Angles (deg)');
save_nii(nii, fullfile(ResFld, ['Sess-' Session '_ang.nii']));
gzip(fullfile(ResFld, ['Sess-' Session '_ang.nii']));
delete(fullfile(ResFld, ['Sess-' Session '_ang.nii']));
% ecc ---
fprintf('Ecc ');
nii = make_nii(result.ecc,[1 1 1],[],[],...
    'pRF fit: Eccentricity (pix)');
save_nii(nii, fullfile(ResFld, ['Sess-' Session '_ecc.nii']));
gzip(fullfile(ResFld, ['Sess-' Session '_ecc.nii']));
delete(fullfile(ResFld, ['Sess-' Session '_ecc.nii']));
% size ---
fprintf('Size ');
nii = make_nii(result.rfsize,[1 1 1],[],[],...
    'pRF fit: RF size (pix)');
save_nii(nii, fullfile(ResFld, ['Sess-' Session '_rfsize.nii']));
gzip(fullfile(ResFld, ['Sess-' Session '_rfsize.nii']));
delete(fullfile(ResFld, ['Sess-' Session '_rfsize.nii']));
% R^2 Goodness of fit ---
fprintf('R2 ');
nii = make_nii(result.R2,[1 1 1],[],[],...
    'pRF fit: R2 Goodnes off fit');
save_nii(nii, fullfile(ResFld, ['Sess-' Session '_R2.nii']));
gzip(fullfile(ResFld, ['Sess-' Session '_R2.nii']));
delete(fullfile(ResFld, ['Sess-' Session '_R2.nii']));

fprintf('>> Done!\n');

