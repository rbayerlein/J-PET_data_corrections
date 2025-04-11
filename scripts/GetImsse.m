function sse = GetImsse(dataA, dataB, c)
%GETIMSSE Summary of this function goes here
%   Detailed explanation goes here
% Get IMage Summed, Squared Error

sse = sum((dataA - c*dataB).^2,'all');

end
