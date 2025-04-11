function [data_out,V] = Mumap2Matmap_JPET(data_in,fname_water, img_size)
%CT2MATMAP Summary of this function goes here
%   Detailed explanation goes here

fprintf('Processing %s...\n',data_in);

disp('Converting u-map to water map...');
disp('Make sure enough worker threads are running for maximum performance!');
t0 = tic; % timer

mumap_in = fread(fopen(data_in, 'rb'), inf, 'float');

data_out = Umap2Watermap(mumap_in,fname_water); % convert u-map to water map
fprintf('Elapsed time is %f seconds.\n',toc(t0));

disp('Converting water map to material index map...');
data_out = Watermap2Matmap(data_out); % convert water map to material index map
% fid=fopen('after_watermap2matmap.img', 'wb');
% fwrite(fid, data_out, 'float');
% fclose(fid);

% next step not needed as input mumap is already in in-house format.
% disp('Flipping image to Eric''s format...');
% data_out = Flip2Eric(data_out); % flip image to Eric's format
% fid=fopen('after_flip2eric.img', 'wb');
% fwrite(fid, data_out, 'float');
% fclose(fid);

data_out = reshape(data_out, img_size);
disp('Flipping image to SimSET''s format...');
data_out = Eric2phg(data_out);
% data_out = rot90(data_out,3);
% fid=fopen('after_rot90.img', 'wb');
% fwrite(fid, data_out, 'float');
% fclose(fid);

end
