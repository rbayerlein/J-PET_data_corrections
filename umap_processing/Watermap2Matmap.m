function data_out = Watermap2Matmap(data_in)
%WATERMAP2MATMAP Summary of this function goes here
%   Detailed explanation goes here

% convert to material index
data_out = data_in + 100; % water00001 starts at 100

% prevent potential bugs since water00001 is the lower water limit
data_out(data_out == 100) = 0; % convert water00001 to air
data_out(data_out > 2301) = 2301; % catch potential out-of-range bug
end
