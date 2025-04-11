function data_out = MaskTxFov(data_in,dcm_info,tx_fov)
%MASKTXFOV Summary of this function goes here
%   Detailed explanation goes here

% this function removes uEXPLORER CT ring artifact
% this function assumes no CT alignment issues

if nargin ~= 3
    error('Usage: MaskTxFov(data_in,dcm_info,tx_fov)');
end

pixel_size = dcm_info.(dicomlookup('0028','0030')); % mm x mm

% see http://matlab.wikia.com/wiki/FAQ#How_do_I_create_a_circle.3F
% and https://www.mathworks.com/matlabcentral/answers/346377-circle-in-a-matrix
% for more info
% assumes image size and pixel size in x/y directions are the same
r = floor(tx_fov/2/pixel_size(1)); % radius (in pixels) - smaller radius preferred to minimize CT ring artifact
x_ctr = ceil(size(data_in,1)/2);
y_ctr = ceil(size(data_in,1)/2);

[n_col,n_row] = meshgrid(1:size(data_in,1) ,1:size(data_in,1));

mask = (n_row - y_ctr).^2 + (n_col - x_ctr).^2 <= r.^2;
mask = repmat(mask,[1,1,size(data_in,3)]);
data_out = data_in .* mask;

end
