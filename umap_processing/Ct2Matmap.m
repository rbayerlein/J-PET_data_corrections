function [data_out,V] = Ct2Matmap(dir_in,hu_limits,fname_water,img_size_out,vox_size_out,tx_fov)
%CT2MATMAP Summary of this function goes here
%   Detailed explanation goes here

if nargin ~= 6
    error('Usage: Ct2Matmap(dir_in,hu_limits,fname_water,img_size_out,vox_size_out,tx_fov)');
end

fprintf('Processing %s...\n',dir_in);

disp('Converting CT to u-map...');
[data_out,V] = Ct2Umap(dir_in,hu_limits); % convert CT to u-map

disp('Applying binary mask...');
dcm_info = dicominfo(sprintf('%s/00000001.dcm',dir_in));
data_out = MaskTxFov(data_out,dcm_info,tx_fov);

disp('Downscaling u-map...');
data_out = DownscaleVolume(data_out,dcm_info,img_size_out,vox_size_out);

disp('Converting u-map to water map...');
disp('Make sure enough worker threads are running for maximum performance!');
t0 = tic; % timer
data_out = Umap2Watermap(data_out,fname_water); % convert u-map to water map
fprintf('Elapsed time is %f seconds.\n',toc(t0));

disp('Converting water map to material index map...');
data_out = Watermap2Matmap(data_out); % convert water map to material index map

%disp('Flipping image to Eric''s format...');
%data_out = Flip2Eric(data_out); % flip image to Eric's format

%disp('Flipping image to SimSET''s format...');
%data_out = Eric2phg(data_out);

end
