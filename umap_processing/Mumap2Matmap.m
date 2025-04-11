function [data_out,V] = Mumap2Matmap(data_in,fname_water)
%CT2MATMAP Summary of this function goes here
%   Detailed explanation goes here

if nargin ~= 2
    error('Usage: Mumap2Matmap(data_in,fname_water)');
end

fprintf('Processing %s...\n',data_in);

disp('Converting u-map to water map...');
disp('Make sure enough worker threads are running for maximum performance!');
t0 = tic; % timer
mumap_in = fread(fopen(data_in, 'rb'), inf, 'float');
%mumap_in = fread(fopen(data_in, 'rb'), inf, 'float');

data_out = Umap2Watermap(mumap_in,fname_water); % convert u-map to water map
fprintf('Elapsed time is %f seconds.\n',toc(t0));

disp('Converting water map to material index map...');
data_out = Watermap2Matmap(data_out); % convert water map to material index map

% next step not needed as input mumap is already in in-house format.
% disp('Flipping image to Eric''s format...');
% data_out = Flip2Eric(data_out); % flip image to Eric's format

disp('Flipping image to SimSET''s format...');
data_out=reshape(data_out,200,200,200);

data_out = Eric2phg(data_out);

end
