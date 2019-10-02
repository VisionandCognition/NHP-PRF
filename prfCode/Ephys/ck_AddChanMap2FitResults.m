subsess = {...
    'Lick','20180807_B2';...
    'Aston','20181004_B1'};

% locations 
base_path = pwd;
if ismac % laptop & portable HDD
    data_path = '/Volumes/MRI_2/PRF_EPHYS';
    home_path = '/Users/chris';
else %
    data_path = '/media/NETDISKS/VCNIN/PRF_EPHYS';
    home_path = '/home/chris';
end
data_fld = fullfile(data_path,'Data_proc');
ch_fld = fullfile(data_path,'Channelmaps');

for su = 1:size(subsess,1)
    fprintf(['=== Adding channelmap to ' subsess{su,1} ...
        ', ' subsess{su,2} ' ===\n']);
    % Load chanmap
    fprintf('Loading channel map\n');
    if strcmp(subsess{su,1},'Lick')
        C=load(fullfile(ch_fld,'channel_area_mapping_lick.mat'));
    elseif strcmp(subsess{su,1},'Aston')
        C=load(fullfile(ch_fld,'channel_area_mapping_aston.mat'));
    end
    
    cd(fullfile(data_fld,subsess{su,1},subsess{su,2},'LFP'));
    fprintf('Appending map to PRF fits\n')
    save([subsess{su,1} '_' subsess{su,2} '_allarrays_PRFLFP'],...
        '-append','C');
    for a=1:8
        fprintf(['Array ' num2str(a) '\n']);
        save([subsess{su,1} '_' subsess{su,2} '_array_' num2str(a) ...
            '_PRFLFP'],'-append','C');
        save([subsess{su,1} '_' subsess{su,2} '_array_' num2str(a) ...
            '_mLFP'],'-append','C');
    end
    cd ../MUA
    fprintf('Appending map to MUA fits\n')
    save([subsess{su,1} '_' subsess{su,2} '_allarrays_PRFMUA'],...
        '-append','C');
    for a=1:8
        fprintf(['Array ' num2str(a) '\n']);
        save([subsess{su,1} '_' subsess{su,2} '_array_' num2str(a) ...
            '_PRFMUA'],'-append','C');
        save([subsess{su,1} '_' subsess{su,2} '_array_' num2str(a) ...
            '_mMUA'],'-append','C');
    end
end
cd(base_path);