function data_out = Eric2phg(data_in)
%ERIC2PHG Summary of this function goes here
%   Detailed explanation goes here

%flip image to SimSET's format
data_out = permute(data_in,[1,2,3]);
data_out = flip(data_out, 2);

end
