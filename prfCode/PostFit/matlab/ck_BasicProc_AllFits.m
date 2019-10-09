function ck_BasicProc_AllFits

% this script will do some basic processing like calculate 
% - means from two-way crossvalidated results
% - X,Y coordinates from polar coordinates
% - highest R2 for crossval
% - distance between crossval runs in location and size

%% load ===================================================================
fitres_path = ...
    '/Users/chris/Dropbox/CURRENT_PROJECTS/NHP_MRI/Projects/pRF/FitResults/';
fprintf('Loading data...');
load(fullfile(fitres_path,'MultiModal','AllFits_cv1'),'R_MRI','R_EPHYS');
fprintf('DONE\n')

%% mri --------------------------------------------------------------------
fprintf('MRI data ---------------------\n');
for r = 1:length(R_MRI) % animals
    for m = 1:length(R_MRI(r).model) % model fits
        fprintf(['Processing MODEL ' R_MRI(r).model(m).prfmodel '\n']);
        
        % convert pixels to degrees
        R_MRI(r).model(m).rfs = R_MRI(r).model(m).rfs./10;
        R_MRI(r).model(m).ecc = R_MRI(r).model(m).ecc./10;
        R_MRI(r).model(m).fwhm = R_MRI(r).model(m).fwhm./10;
        
        F={'rfs','gain','R2','fwhm','sdratio','normamp','expt'};
        
        % avg some values
        for f=1:length(F)
            if ~isempty(R_MRI(r).model(m).(F{f}))
                R_MRI(r).model(m).avg.(F{f}) = mean(R_MRI(r).model(m).(F{f}),1);
            end
        end
        
        % XY
        R_MRI(r).model(m).X = R_MRI(r).model(m).ecc.*cosd(R_MRI(r).model(m).ang);
        R_MRI(r).model(m).Y = R_MRI(r).model(m).ecc.*sind(R_MRI(r).model(m).ang);
        
        % avg XY
        R_MRI(r).model(m).avg.X=mean(R_MRI(r).model(m).X,1);
        R_MRI(r).model(m).avg.Y=mean(R_MRI(r).model(m).Y,1);
        R_MRI(r).model(m).avg.ecc=sqrt( ...
            R_MRI(r).model(m).avg.X.^2 + R_MRI(r).model(m).avg.Y.^2);
        
        % avg angle
        rawang = atand(R_MRI(r).model(m).avg.Y./R_MRI(r).model(m).avg.X);
        rawang(R_MRI(r).model(m).avg.X<0 & R_MRI(r).model(m).avg.Y<0) = ...
            rawang(R_MRI(r).model(m).avg.X<0 & R_MRI(r).model(m).avg.Y<0)-180;
        rawang(R_MRI(r).model(m).avg.X<0 & R_MRI(r).model(m).avg.Y>=0) = ...
            rawang(R_MRI(r).model(m).avg.X<0 & R_MRI(r).model(m).avg.Y>=0)+180;
        R_MRI(r).model(m).avg.ang = rawang;
        
        % dXY
        R_MRI(r).model(m).diff.XY = sqrt(...
            diff(R_MRI(r).model(m).X,1,2).^2 + ...
            diff(R_MRI(r).model(m).Y,1,2).^2);
        
        % dRFS
        R_MRI(r).model(m).diff.rfs = diff(R_MRI(r).model(m).rfs,1,2);
        
        % get maximum R2 version of some values
        [R_MRI(r).model(m).max.R2, R_MRI(r).model(m).max.R2_idx] = ...
            max(R_MRI(r).model(m).R2,[],1);
        R_MRI(r).model(m).max.R2_idx = ...
            [R_MRI(r).model(m).max.R2_idx==1; R_MRI(r).model(m).max.R2_idx==2];
        
        F={'rfs','gain','ecc','ang','fwhm','sdratio','normamp','expt'};
        for f=1:length(F)
            if ~isempty(R_MRI(r).model(m).(F{f}))
                R_MRI(r).model(m).max.(F{f}) = ...
                    R_MRI(r).model(m).(F{f})(R_MRI(r).model(m).max.R2_idx)';
            end
        end
        R_MRI(r).model(m).max.X = R_MRI(r).model(m).X(R_MRI(r).model(m).max.R2_idx)';
        R_MRI(r).model(m).max.Y = R_MRI(r).model(m).Y(R_MRI(r).model(m).max.R2_idx)';
    end
end

%% ephys ------------------------------------------------------------------
fprintf('Ephys data ---------------------\n');
for r = 1:length(R_EPHYS) % animals
    for m = 1:length(R_EPHYS(r).model) % model fits
        fprintf(['Processing MODEL ' R_EPHYS(r).model(m).prfmodel '\n']);
        if ~strcmp(R_EPHYS(r).model(m).prfmodel,'classicRF')
            
            %% MUA --------------------------------------------------------
            fprintf('MUA ');
            for i = 1:8 % all instances
                % make name consistent with MRI
                R_EPHYS(r).model(m).MUA(i).rfs = R_EPHYS(r).model(m).MUA(i).rfsize;
                
                % convert pixels to degrees
                R_EPHYS(r).model(m).MUA(i).ecc = ...
                    R_EPHYS(r).model(m).MUA(i).ecc .* R_EPHYS(r).model(m).MUA(i).Pix2Deg;
                R_EPHYS(r).model(m).MUA(i).rfs = ...
                    R_EPHYS(r).model(m).MUA(i).rfs .* R_EPHYS(r).model(m).MUA(i).Pix2Deg;
                
                % calculate fwhm
                R_EPHYS(r).model(m).MUA(i).fwhm = ...
                    R_EPHYS(r).model(m).MUA(i).rfs.*(2*sqrt(2*log(2)));
                              
                F={'rfs','gain','R2','fwhm','sdratio','normamp','expt'};
                
                for f=1:length(F)
                    if isfield(R_EPHYS(r).model(m).MUA(i), F{f})
                        R_EPHYS(r).model(m).MUA(i).avg.(F{f}) = mean(R_EPHYS(r).model(m).MUA(i).(F{f}),2);
 
                    end
                end
                
                % XY
                R_EPHYS(r).model(m).MUA(i).X = ...
                    R_EPHYS(r).model(m).MUA(i).ecc.*cosd(R_EPHYS(r).model(m).MUA(i).ang);
                R_EPHYS(r).model(m).MUA(i).Y = ...
                    R_EPHYS(r).model(m).MUA(i).ecc.*sind(R_EPHYS(r).model(m).MUA(i).ang);
                
                % avg XY
                R_EPHYS(r).model(m).MUA(i).avg.X=mean(R_EPHYS(r).model(m).MUA(i).X,2)';
                R_EPHYS(r).model(m).MUA(i).avg.Y=mean(R_EPHYS(r).model(m).MUA(i).Y,2)';
                R_EPHYS(r).model(m).MUA(i).avg.ecc=sqrt( ...
                    R_EPHYS(r).model(m).MUA(i).avg.X.^2 + R_EPHYS(r).model(m).MUA(i).avg.Y.^2);
                
                % avg angle
                rawang = atand(R_EPHYS(r).model(m).MUA(i).avg.Y./R_EPHYS(r).model(m).MUA(i).avg.X);
                rawang(R_EPHYS(r).model(m).MUA(i).avg.X<0 & R_EPHYS(r).model(m).MUA(i).avg.Y<0) = ...
                    rawang(R_EPHYS(r).model(m).MUA(i).avg.X<0 & R_EPHYS(r).model(m).MUA(i).avg.Y<0)-180;
                rawang(R_EPHYS(r).model(m).MUA(i).avg.X<0 & R_EPHYS(r).model(m).MUA(i).avg.Y>=0) = ...
                    rawang(R_EPHYS(r).model(m).MUA(i).avg.X<0 & R_EPHYS(r).model(m).MUA(i).avg.Y>=0)+180;
                R_EPHYS(r).model(m).MUA(i).avg.ang = rawang;
                
                % dXY
                R_EPHYS(r).model(m).MUA(i).diff.XY = sqrt(...
                    diff(R_EPHYS(r).model(m).MUA(i).X,1,2).^2 + ...
                    diff(R_EPHYS(r).model(m).MUA(i).Y,1,2).^2);
                % dRFS
                R_EPHYS(r).model(m).MUA(i).diff.rfs = diff(R_EPHYS(r).model(m).MUA(i).rfs,1,2);
                
                % get maximum R2 version of some values
                [R_EPHYS(r).model(m).MUA(i).max.R2, R_EPHYS(r).model(m).MUA(i).max.R2_idx] = ...
                    max(R_EPHYS(r).model(m).MUA(i).R2,[],2);
                R_EPHYS(r).model(m).MUA(i).max.R2_idx = ...
                     [R_EPHYS(r).model(m).MUA(i).max.R2_idx==1;...
                      R_EPHYS(r).model(m).MUA(i).max.R2_idx==2];
                
                F={'rfs','gain','ecc','ang','fwhm','sdratio','normamp','expt'};
                for f=1:length(F)
                    if isfield(R_EPHYS(r).model(m).MUA(i), F{f}) && ~isempty(R_EPHYS(r).model(m).MUA(i).(F{f}))
                        R_EPHYS(r).model(m).MUA(i).max.(F{f}) = ...
                            R_EPHYS(r).model(m).MUA(i).(F{f})(R_EPHYS(r).model(m).MUA(i).max.R2_idx)';
                    end
                end
                R_EPHYS(r).model(m).MUA(i).max.X = R_EPHYS(r).model(m).MUA(i).X(R_EPHYS(r).model(m).MUA(i).max.R2_idx)';
                R_EPHYS(r).model(m).MUA(i).max.Y = R_EPHYS(r).model(m).MUA(i).Y(R_EPHYS(r).model(m).MUA(i).max.R2_idx)';

            end
            if isfield(R_EPHYS(r).model(m).MUA,'rfsize')
                R_EPHYS(r).model(m).MUA = rmfield(R_EPHYS(r).model(m).MUA,'rfsize');
            end
            
            %% LFP --------------------------------------------------------
            fprintf('LFP\n');
            for fb=1:5
                for i = 1:8 % all instances
                    % make name consistent with MRI
                    R_EPHYS(r).model(m).LFP(i,fb).rfs = R_EPHYS(r).model(m).LFP(i,fb).rfsize;
                    
                    % convert pixels to degrees
                    R_EPHYS(r).model(m).LFP(i,fb).ecc = ...
                        R_EPHYS(r).model(m).LFP(i,fb).ecc .* R_EPHYS(r).model(m).LFP(i,fb).Pix2Deg;
                    R_EPHYS(r).model(m).LFP(i,fb).rfs = ...
                        R_EPHYS(r).model(m).LFP(i,fb).rfs .* R_EPHYS(r).model(m).LFP(i,fb).Pix2Deg;
                    
                    % calculate fwhm
                    R_EPHYS(r).model(m).LFP(i,fb).fwhm = ...
                        R_EPHYS(r).model(m).LFP(i,fb).rfs.*(2*sqrt(2*log(2)));
                
                    F={'rfs','gain','R2','fwhm','sdratio','normamp','expt'};
                    
                    for f=1:length(F)
                        if isfield(R_EPHYS(r).model(m).LFP(i,fb), F{f})
                            R_EPHYS(r).model(m).LFP(i,fb).avg.(F{f}) = mean(R_EPHYS(r).model(m).LFP(i,fb).(F{f}),2)';
                        end
                    end
                    
                    % XY
                    R_EPHYS(r).model(m).LFP(i,fb).X = ...
                        R_EPHYS(r).model(m).LFP(i,fb).ecc.*cosd(R_EPHYS(r).model(m).LFP(i,fb).ang);
                    R_EPHYS(r).model(m).LFP(i,fb).Y = ...
                        R_EPHYS(r).model(m).LFP(i,fb).ecc.*sind(R_EPHYS(r).model(m).LFP(i,fb).ang);
                    
                    % avg XY
                    R_EPHYS(r).model(m).LFP(i,fb).avg.X=mean(R_EPHYS(r).model(m).LFP(i,fb).X,2)';
                    R_EPHYS(r).model(m).LFP(i,fb).avg.Y=mean(R_EPHYS(r).model(m).LFP(i,fb).Y,2)';
                    R_EPHYS(r).model(m).LFP(i,fb).avg.ecc=sqrt( ...
                        R_EPHYS(r).model(m).LFP(i,fb).avg.X.^2 + R_EPHYS(r).model(m).LFP(i,fb).avg.Y.^2);
                    
                    % avg angle
                    rawang = atand(R_EPHYS(r).model(m).LFP(i,fb).avg.Y./R_EPHYS(r).model(m).LFP(i,fb).avg.X);
                    rawang(R_EPHYS(r).model(m).LFP(i,fb).avg.X<0 & R_EPHYS(r).model(m).LFP(i,fb).avg.Y<0) = ...
                        rawang(R_EPHYS(r).model(m).LFP(i,fb).avg.X<0 & R_EPHYS(r).model(m).LFP(i,fb).avg.Y<0)-180;
                    rawang(R_EPHYS(r).model(m).LFP(i,fb).avg.X<0 & R_EPHYS(r).model(m).LFP(i,fb).avg.Y>=0) = ...
                        rawang(R_EPHYS(r).model(m).LFP(i,fb).avg.X<0 & R_EPHYS(r).model(m).LFP(i,fb).avg.Y>=0)+180;
                    R_EPHYS(r).model(m).LFP(i,fb).avg.ang = rawang;
                    
                    % dXY
                    R_EPHYS(r).model(m).LFP(i,fb).diff.XY = sqrt(...
                        diff(R_EPHYS(r).model(m).LFP(i,fb).X,1,2).^2 + ...
                        diff(R_EPHYS(r).model(m).LFP(i,fb).Y,1,2).^2);
                    % dRFS
                    R_EPHYS(r).model(m).LFP(i,fb).diff.rfs = diff(R_EPHYS(r).model(m).LFP(i,fb).rfs,1,2);
                    
                    % get maximum R2 version of some values
                    [R_EPHYS(r).model(m).LFP(i,fb).max.R2, R_EPHYS(r).model(m).LFP(i,fb).max.R2_idx] = ...
                        max(R_EPHYS(r).model(m).LFP(i,fb).R2,[],2);
                    R_EPHYS(r).model(m).LFP(i,fb).max.R2_idx = ...
                        [R_EPHYS(r).model(m).LFP(i,fb).max.R2_idx==1;...
                        R_EPHYS(r).model(m).LFP(i,fb).max.R2_idx==2];
                    
                    F={'rfs','gain','ecc','ang','fwhm','sdratio','normamp','expt'};
                    for f=1:length(F)
                        if isfield(R_EPHYS(r).model(m).LFP(i,fb), F{f}) && ~isempty(R_EPHYS(r).model(m).LFP(i,fb).(F{f}))
                            R_EPHYS(r).model(m).LFP(i,fb).max.(F{f}) = ...
                                R_EPHYS(r).model(m).LFP(i,fb).(F{f})(R_EPHYS(r).model(m).LFP(i,fb).max.R2_idx)';
                        end
                    end
                    R_EPHYS(r).model(m).LFP(i,fb).max.X = R_EPHYS(r).model(m).LFP(i,fb).X(R_EPHYS(r).model(m).LFP(i,fb).max.R2_idx)';
                    R_EPHYS(r).model(m).LFP(i,fb).max.Y = R_EPHYS(r).model(m).LFP(i,fb).Y(R_EPHYS(r).model(m).LFP(i,fb).max.R2_idx)';
                end
            end
            if isfield(R_EPHYS(r).model(m).LFP,'rfsize')
                R_EPHYS(r).model(m).LFP = rmfield(R_EPHYS(r).model(m).LFP,'rfsize');
            end
            
        else % classic RF mapping
            fprintf('MUA\n');
         
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
        end
    end 
end