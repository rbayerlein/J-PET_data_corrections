function data_out = Umap2Watermap(data_in,fname_water)
%UMAP2WATERMAP Summary of this function goes here
%   Detailed explanation goes here

water_in = fread(fopen(fname_water),'double');
water_in = water_in/10;
data_in = data_in/10; % convert from cm^-1 to mm^-1

data_out = zeros(size(data_in)); % initialization


% convert to water density
parfor i = 1:numel(data_in)
    [~,data_out(i)] = min(abs(data_in(i)-water_in)); % save index of u0
end
data_out = floor(data_out/10);

end