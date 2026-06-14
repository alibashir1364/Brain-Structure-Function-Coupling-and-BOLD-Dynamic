% --------------------------------------------------------
% V : 3D volume of Glasser atlas
%     voxel values:
%       0   -> background
%       1-360 -> Glasser ROI labels
% --------------------------------------------------------
V = niftiread('...\data\volume\volume\HCP\100307\HCPMMP1_standard.nii');

numROIs = 360;

roi_volume = zeros(numROIs, 1);

% ----------------------------------------

% ----------------------------------------
for roi = 1:numROIs
    
   
    roi_volume(roi) = nnz(V == roi);
    
end

% ----------------------------------------

% HCP volume: 2mm isotropic
% ----------------------------------------
voxel_volume_mm3 = 2 * 2 * 2;

roi_volume_mm3 = roi_volume * voxel_volume_mm3;

