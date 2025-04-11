function data_out = SaveImg2Amap_original(fname_in,fname_amap,fname_phg_act_table,n_bins,img_size,vox_size)
%SAVEIMG2AMAP Summary of this function goes here
%   Detailed explanation goes here

if nargin ~= 6
    error('Usage: SaveImg2Amap(fname_in,fname_amap,fname_phg_act_table,n_bins,img_size,vox_size)');
end
pause(0.1);


fprintf('Processing %s...\n',fname_in);
fprintf('Image size: %d, %d, %d', img_size(1), img_size(2), img_size(3));

fid_data_in = fopen(fname_in, 'rb');
pause(1);

data_in = fread(fid_data_in, 'float'); % a.u.
pause(1);

s=size(data_in);
data_in = reshape(data_in,img_size);
size(data_in)

% scale entire activity image to 10 mCi equivalent ("QF")
% quick tip: 1 mCi in 2m pipe is about 1 kBq/cc
% voxel activity concentration (Bq/ml) * voxel size (ml) = voxel activity (Bq)
% image activity (Bq) can be rewritten as sum(data_in)*prod(vox_size)/1000
% prod(vox_size)/1000 converts mm3 to cm3
% sum(data_in)*prod(vox_size)/1000 (Bq) * QF (Ci/Bq) = 0.01 Ci (i.e. 10 mCi)

qf = 0.01/(sum(data_in,'all')*prod(vox_size)/1000);

% apply QF to get activity concentration
data_out = data_in*qf; % Ci/cc

% activity concentration resolution
ac_res = max(data_out,[],'all')/n_bins;

% convert to amap indecies
data_out = floor(data_out/ac_res); % prevent out-of-bound indecies
data_out(data_out == n_bins) = n_bins-1; % prevent out-of-bound indecies

%% write out activity map

% set +/- Inf and NaN values to 0 to prevent errors
data_out(isinf(data_out)) = 0;
data_out(isnan(data_out)) = 0;

% flip to SimSET's format
data_out = Eric2phg(data_out);
fid_amap = fopen(fname_amap,'w');
fwrite(fid_amap, data_out,'int32');

%% write out phg_act_table

fid_phg_act_table = fopen(fname_phg_act_table,'w');
fprintf(fid_phg_act_table,'%d\n',n_bins);
fprintf(fid_phg_act_table,'%e\n',ac_res*(0:n_bins-1));
fclose(fid_amap);
fclose(fid_phg_act_table);
end
