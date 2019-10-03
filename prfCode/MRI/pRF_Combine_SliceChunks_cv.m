function pRF_Combine_SliceChunks_cv(Monkey,Session)
% For prf analysis that is run on LISA in separate slice chunks, this
% function recombines the results in volumes

tool_basepath = '/Users/chris/Dropbox/MATLAB_NONGIT/TOOLBOX';
addpath(genpath(fullfile(tool_basepath, 'NIfTI')));

homefld = pwd;
ResFld =  ['~/Documents/MRI_ANALYSIS/NHP-PRF/FitResults/MRI/' Monkey '/' Session];
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
        if isfield(result,'expt')
            result.expt     = cat(3,result.expt,R.result.expt);
        end
        if isfield(result,'sdratio')
            result.sdratio     = cat(3,result.sdratio,R.result.sdratio);
        end
        if isfield(result,'normamp')
            result.normamp     = cat(3,result.normamp,R.result.normamp);
        end
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

% fwhm
fprintf('Calculating FWHM..\n');
if isfield(result,'sdratio')
    result.fwhm=zeros(size(result.rfsize));
    for x=1:size(result.rfsize,1)
        for y=1:size(result.rfsize,2)
            for z=1:size(result.rfsize,3)
                for cv = 1:size(result.rfsize,4)
                    sig1 = result.rfsize(x,y,z,cv);
                    sig2 = result.sdratio(x,y,z,cv).*result.rfsize(x,y,z,cv);
                    a2 = result.normamp(x,y,z,cv);
                    p=-100:0.1:100;
                    g = normpdf(p,0,sig1) - a2*normpdf(p,0,sig2);
                    pp = p(g>=(max(g)/2));
                    if ~isempty(pp)
                        result.fwhm(x,y,z,cv) = pp(end)-pp(1);
                    else
                        result.fwhm(x,y,z,cv) = 0;
                    end
                end
            end
        end
    end
else
    for cv = 1:size(result.rfsize,4)
        result.fwhm(:,:,:,cv) = result.rfsize(:,:,:,cv).*(2*sqrt(2*log(2)));
    end
end

% save the result ----
fprintf('Merging and saving the result: mat-file, ');
save(fullfile(ResFld,['pRF_Sess-' Session]),'result','-v7.3');

% also save as nifti files
for cv=1:2
    % angle ---
    fprintf('Angles, ');
    nii = make_nii(result.ang(:,:,:,cv),[1 1 1],[],16,...
        'pRF fit: Angles (deg)');
    nii.hdr.hist = result.hdr_ref.hist;
    nii.hdr.dime.datatype = 64; nii.hdr.dime.bitpix = 64;
    save_nii(nii, fullfile(ResFld, ...
        ['Sess-' Session '_ang_' num2str(cv) '.nii']));
    gzip(fullfile(ResFld, ['Sess-' Session '_ang_' num2str(cv) '.nii']));
    delete(fullfile(ResFld, ['Sess-' Session '_ang_' num2str(cv) '.nii']));
    
    % split to real(-cos)/imag(sin) components
    nii = make_nii(-cosd(result.ang(:,:,:,cv)),[1 1 1],[],16,...
        'pRF fit: Angles_real (deg)');
    nii.hdr.hist = result.hdr_ref.hist;
    nii.hdr.dime.datatype = 64; nii.hdr.dime.bitpix = 64;
    save_nii(nii, fullfile(ResFld, ...
        ['Sess-' Session '_real_' num2str(cv) '.nii']));
    gzip(fullfile(ResFld, ['Sess-' Session '_real_' num2str(cv) '.nii']));
    delete(fullfile(ResFld, ['Sess-' Session '_real_' num2str(cv) '.nii']));
    nii = make_nii(sind(result.ang(:,:,:,cv)),[1 1 1],[],16,...
        'pRF fit: Angles_imag (deg)');
    nii.hdr.hist = result.hdr_ref.hist;
    nii.hdr.dime.datatype = 64; nii.hdr.dime.bitpix = 64;
    save_nii(nii, fullfile(ResFld, ...
        ['Sess-' Session '_imag_' num2str(cv) '.nii']));
    gzip(fullfile(ResFld, ['Sess-' Session '_imag_' num2str(cv) '.nii']));
    delete(fullfile(ResFld, ['Sess-' Session '_imag_' num2str(cv) '.nii']));
    
    % ecc ---
    fprintf('Ecc, ');
    nii = make_nii(result.ecc(:,:,:,cv),[1 1 1],[],16,...
        'pRF fit: Eccentricity (pix)');
    nii.hdr.hist = result.hdr_ref.hist;
    nii.hdr.dime.datatype = 64; nii.hdr.dime.bitpix = 64;
    save_nii(nii, fullfile(ResFld, ...
        ['Sess-' Session '_ecc_' num2str(cv) '.nii']));
    gzip(fullfile(ResFld, ['Sess-' Session '_ecc_' num2str(cv) '.nii']));
    delete(fullfile(ResFld, ['Sess-' Session '_ecc_' num2str(cv) '.nii']));
    
    % size ---
    fprintf('Size, ');
    nii = make_nii(result.rfsize(:,:,:,cv),[1 1 1],[],16,...
        'pRF fit: RF size (pix)');
    nii.hdr.hist = result.hdr_ref.hist;
    nii.hdr.dime.datatype = 64; nii.hdr.dime.bitpix = 64;
    save_nii(nii, fullfile(ResFld, ...
        ['Sess-' Session '_rfsize_' num2str(cv) '.nii']));
    gzip(fullfile(ResFld, ['Sess-' Session '_rfsize_' num2str(cv) '.nii']));
    delete(fullfile(ResFld, ['Sess-' Session '_rfsize_' num2str(cv) '.nii']));
    
    % exponential ---
    if isfield(result,'expt')
        fprintf('Expt, ');
        nii = make_nii(result.expt(:,:,:,cv),[1 1 1],[],16,...
            'pRF fit: RF size (pix)');
        nii.hdr.hist = result.hdr_ref.hist;
        nii.hdr.dime.datatype = 64; nii.hdr.dime.bitpix = 64;
        save_nii(nii, fullfile(ResFld, ...
            ['Sess-' Session '_expt_' num2str(cv) '.nii']));
        gzip(fullfile(ResFld, ['Sess-' Session '_expt_' num2str(cv) '.nii']));
        delete(fullfile(ResFld, ['Sess-' Session '_expt_' num2str(cv) '.nii']));
        nii = make_nii(result.expt(:,:,:,cv).^2,[1 1 1],[],16,...
            'pRF fit: RF size (pix)');
        nii.hdr.hist = result.hdr_ref.hist;
        nii.hdr.dime.datatype = 64; nii.hdr.dime.bitpix = 64;
        save_nii(nii, fullfile(ResFld, ...
            ['Sess-' Session '_exptsq_' num2str(cv) '.nii']));
        gzip(fullfile(ResFld, ['Sess-' Session '_exptsq_' num2str(cv) '.nii']));
        delete(fullfile(ResFld, ['Sess-' Session '_exptsq_' num2str(cv) '.nii']));
    end
    
    % sdratio
    if isfield(result,'sdratio')
        fprintf('sdratio, ');
        nii = make_nii(result.sdratio(:,:,:,cv),[1 1 1],[],16,...
            'pRF fit: sdratio');
        nii.hdr.hist = result.hdr_ref.hist;
        nii.hdr.dime.datatype = 64; nii.hdr.dime.bitpix = 64;
        save_nii(nii, fullfile(ResFld, ...
            ['Sess-' Session '_sdratio_' num2str(cv) '.nii']));
        gzip(fullfile(ResFld, ['Sess-' Session '_sdratio_' num2str(cv) '.nii']));
        delete(fullfile(ResFld, ['Sess-' Session '_sdratio_' num2str(cv) '.nii']));
        fprintf('sd2, ');
        sd2 = result.sdratio(:,:,:,cv).*result.rfsize(:,:,:,cv);
        nii = make_nii(sd2,[1 1 1],[],16,...
            'pRF fit: sd2');
        nii.hdr.hist = result.hdr_ref.hist;
        nii.hdr.dime.datatype = 64; nii.hdr.dime.bitpix = 64;
        save_nii(nii, fullfile(ResFld, ...
            ['Sess-' Session '_sd2_' num2str(cv) '.nii']));
        gzip(fullfile(ResFld, ['Sess-' Session '_sd2_' num2str(cv) '.nii']));
        delete(fullfile(ResFld, ['Sess-' Session '_sd2_' num2str(cv) '.nii']));
        
        % fwhm
        fprintf('FWHM, ');
        nii = make_nii(result.fwhm(:,:,:,cv),[1 1 1],[],16,...
            'pRF fit: FWHM (pix)');
        nii.hdr.hist = result.hdr_ref.hist;
        nii.hdr.dime.datatype = 64; nii.hdr.dime.bitpix = 64;
        save_nii(nii, fullfile(ResFld, ...
            ['Sess-' Session '_FWHM_' num2str(cv) '.nii']));
        gzip(fullfile(ResFld, ['Sess-' Session '_FWHM_' num2str(cv) '.nii']));
        delete(fullfile(ResFld, ['Sess-' Session '_FWHM_' num2str(cv) '.nii']));
    else    
        % fwhm
        fprintf('FWHM, ');
        nii = make_nii(result.fwhm(:,:,:,cv),[1 1 1],[],16,...
            'pRF fit: FWHM (pix)');
        nii.hdr.hist = result.hdr_ref.hist;
        nii.hdr.dime.datatype = 64; nii.hdr.dime.bitpix = 64;
        save_nii(nii, fullfile(ResFld, ...
            ['Sess-' Session '_FWHM_' num2str(cv) '.nii']));
        gzip(fullfile(ResFld, ['Sess-' Session '_FWHM_' num2str(cv) '.nii']));
        delete(fullfile(ResFld, ['Sess-' Session '_FWHM_' num2str(cv) '.nii']));
    end
    
    % normamp
    if isfield(result,'normamp')
        fprintf('normamp, ');
        nii = make_nii(result.normamp(:,:,:,cv),[1 1 1],[],16,...
            'pRF fit: normamp');
        nii.hdr.hist = result.hdr_ref.hist;
        nii.hdr.dime.datatype = 64; nii.hdr.dime.bitpix = 64;
        save_nii(nii, fullfile(ResFld, ...
            ['Sess-' Session '_normamp_' num2str(cv) '.nii']));
        gzip(fullfile(ResFld, ['Sess-' Session '_normamp_' num2str(cv) '.nii']));
        delete(fullfile(ResFld, ['Sess-' Session '_normamp_' num2str(cv) '.nii']));
    end
    
    % R^2 Goodness of fit ---
    fprintf('R2 ');
    nii = make_nii(result.R2(:,:,:,cv),[1 1 1],[],16,...
        'pRF fit: R2 Goodnes off fit');
    nii.hdr.hist = result.hdr_ref.hist;
    nii.hdr.dime.datatype = 64; nii.hdr.dime.bitpix = 64;
    save_nii(nii, fullfile(ResFld, ...
        ['Sess-' Session '_R2_' num2str(cv) '.nii']));
    gzip(fullfile(ResFld, ['Sess-' Session '_R2_' num2str(cv) '.nii']));
    delete(fullfile(ResFld, ['Sess-' Session '_R2_' num2str(cv) '.nii']));
end

fprintf('>> Done!\n');
cd(homefld)
