function [data_out,V] = Ct2Umap_only(dir_in,hu_limits,img_size_out,vox_size_out,tx_fov)
%CT2UMAP_ONLY Summary of this function goes here
%   Detailed explanation goes here

if nargin ~= 5
    error('Usage: Ct2Umap_only(dir_in,hu_limits,img_size_out,vox_size_out,tx_fov)');
end

fprintf('Processing %s...\n',dir_in);

disp('Converting CT to u-map...');
[data_out,V] = Ct2Umap(dir_in,hu_limits); % convert CT to u-map

disp('Applying binary mask...');
dcm_info = dicominfo(sprintf('%s/00000001.dcm',dir_in));
data_out = MaskTxFov(data_out,dcm_info,tx_fov);

disp('Downscaling u-map...');
data_out = DownscaleVolume(data_out,dcm_info,img_size_out,vox_size_out);

disp('Flipping image to Eric''s format...');
data_out = Flip2Eric(data_out); % flip image to Eric's format

end
