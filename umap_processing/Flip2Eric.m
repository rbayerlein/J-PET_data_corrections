function data_out = Flip2Eric(data_in)
%FLIP2ERIC Summary of this function goes here
%   Detailed explanation goes here

%flip image to Eric's format
data_out = flip(data_in,3); % head first
data_out = permute(data_out,[2,1,3]);
data_out = rot90(data_out,2);

end