function [data_out,V] = Ct2Umap(dir_in,hu_limits)
% hu_limits is defined as [hu_min,hu_max]
% hu_min is the threshold value at which any values smaller than hu_min
%   is set to -1000 (air)
% hu_max is the threshold value at which any values greater than hu_max
%   is set to hu_max

% fix dir_in path if needed
if dir_in(end) ~= '/'
    dir_in = strcat(dir_in,'/');
end

if numel(hu_limits) ~= 2
    error('Usage: [data_out,V] = Ct2Umap(dir_in,[hu_min,hu_max])');
end

if hu_limits(1) < -1000
    error('hu_min cannot be smaller than -1000!');
end

% get DICOM info
dcm_info = dicominfo(strcat(dir_in,'00000001.dcm'));

% read CT images
V = double(dicomreadVolume(dir_in)); % dicom default is uint16
V = squeeze(V);

% rescale to HU units
V = V * dcm_info.RescaleSlope + dcm_info.RescaleIntercept;

% limit the minimum HU because attenuation coefficients cannot be negative
V(V < hu_limits(1)) = -1000;

% limit the maximum HU to minimize CT artifacts
V(V > hu_limits(2)) = hu_limits(2);

% make a copy so as to not modify the original
data_out = V;

% convert to u-map
% see Method for transforming CT images for attenuation correction in PET/CT imaging
% for more info.

% check DICOM kVp
% a and b are in units of mm^-1
% bp is the break point (HU)
fprintf('Image kVp is: %d\n',dcm_info.KVP);
if dcm_info.KVP == 140
    a = 5.64e-6;
    b = 4.08e-3;
    bp = 30;
else
    error('Invalid kVp setting!');
end

% below breakpoint
data_out(V < bp) = 9.6e-6 * (V(V < bp) + 1000); % mm^-1

% above breakpoint
data_out(V >= bp) = a * (V(V >= bp) + 1000) + b;  % mm^-1

end