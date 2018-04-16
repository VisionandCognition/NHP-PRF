%imagesc(options.vxs(:,:,30))
voxcor = [68,18,30];
voxel={};
for i=1:size(fmri_data,2)
    voxel{i}=reshape(...
        fmri_data{i}(voxcor(1),voxcor(2),voxcor(3),:),...
        [1 size(fmri_data{i},4)]); % some visual voxel
end
result = analyzePRF(stimulus,voxel,TR);