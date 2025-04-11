function data_out = DownscaleVolume(data_in,dcm_info,img_size_out,vox_size_out)
%DOWNSCALEVOLUME Summary of this function goes here
%   Detailed explanation goes here

if nargin ~= 4
    error('Usage: DownscaleVolume(data_in,dcm_info,img_size_out,vox_size_out)');
end

% this function does not add zero padding
% check image dim to make sure output image is smaller than input image

pixel_size = dcm_info.(dicomlookup('0028','0030')); % mm x mm
vox_size_in = [pixel_size(1),pixel_size(2),dcm_info.(dicomlookup('0018','0050'))]; % mm x mm x mm

xdim_in = size(data_in,1)*vox_size_in(1); % mm
ydim_in = size(data_in,2)*vox_size_in(2); % mm
zdim_in = size(data_in,3)*vox_size_in(3); % mm

xdim_out = img_size_out(1)*vox_size_out(1); % mm
ydim_out = img_size_out(2)*vox_size_out(2); % mm
zdim_out = img_size_out(3)*vox_size_out(3); % mm

if xdim_in < xdim_out || ydim_in < ydim_out
    error('Input image dim (x/y) is smaller than output image dim!');
end

while zdim_in < zdim_out
    disp('WARNING: input image dim (z) is smaller than output image dim!');
    disp('Zero-padding the input images...');
    padsize = 1; % pad padsize*2 slices at a time
    data_in = padarray(data_in,[0,0,padsize],0);
    zdim_in = size(data_in,3)*vox_size_in(3); % mm
end

% create input and query (output) grids (see interp3 documentation for more info)

[X,Y,Z] = meshgrid(-xdim_in/2+vox_size_in(1)/2:vox_size_in(1):xdim_in/2-vox_size_in(1)/2, ...
    -ydim_in/2+vox_size_in(2)/2:vox_size_in(2):ydim_in/2-vox_size_in(2)/2, ...
    -zdim_in/2+vox_size_in(3)/2:vox_size_in(3):zdim_in/2-vox_size_in(3)/2);

[Xq,Yq,Zq] = meshgrid(-xdim_out/2+vox_size_out(1)/2:vox_size_out(1):xdim_out/2 -vox_size_out(1)/2, ...
    -ydim_out/2+vox_size_out(2)/2:vox_size_out(2):ydim_out/2-vox_size_out(2)/2, ...
    -zdim_out/2+vox_size_out(3)/2:vox_size_out(3):zdim_out/2-vox_size_out(3)/2);

data_out = interp3(X,Y,Z,data_in,Xq,Yq,Zq);

end
