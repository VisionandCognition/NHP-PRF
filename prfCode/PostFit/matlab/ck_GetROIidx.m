function roiidx = ck_GetROIidx(roilabels,rois)
% roilabels is a cell of desired ROIs to include    
% rois is a 2-column vector with {label , [idx(s)]}
    
    roiidx = []; roilab = rois(:,1);
    
    for idx = 1:length(roilabels)
        roiidx = cat(2,roiidx,...
            rois{strcmp(roilab, roilabels{idx}),2});
    end
end