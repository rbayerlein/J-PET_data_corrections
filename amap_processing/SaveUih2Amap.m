function data_out = SaveUih2Amap(dir_in,num_frames,frame_num,fname_amap,fname_phg_act_table,n_bins)
%SAVEUIH2AMAP Summary of this function goes here
%   Detailed explanation goes here

if nargin ~= 6
    error('Usage: SaveUih2Amap(fname_in,num_frames,frame_num,fname_amap,fname_phg_act_table,n_bins)');
end

fprintf('Processing %s...\n',dir_in);
disp('Note: activity image is in UIH format!');

% determine file type (2-D vs 3-D)
ls_dir_in = dir(strcat(dir_in,'/','*.dcm')); % equivalent to ls *.dcm
num_files = length(ls_dir_in); % number of files in dir_in
num_files_per_frame = num_files/num_frames;

is_3d_dicom = false; % initialization
if num_files_per_frame == 1
    disp('3-D DICOM files detected.');
    is_3d_dicom = true;
else
    disp('2-D DICOM files detected.');
end

% read dicom info
if is_3d_dicom
    dcm_info = dicominfo(strcat(dir_in,'/',sprintf('%08d',frame_num),'.dcm'));
else
    dcm_info = dicominfo(strcat(dir_in,'/','00000001.dcm'));
end

% stdout
fprintf('Processing UIH frame #%d...\n',frame_num);

% read data
if is_3d_dicom
    data_in = dicomread(dcm_info);
    data_in = squeeze(data_in);
else
    data_in = zeros(dcm_info.Rows,dcm_info.Columns,dcm_info.NumberOfSlices); % initialization
    for i = 1:dcm_info.NumberOfSlices
        fname_in = strcat(dir_in,'/',sprintf('%08d',i),'.dcm');
        data_in(:,:,i) = dicomread(fname_in);
    end
end

% rescale data
data_in = double(data_in); % convert from uint16 to double prior to operations

if is_3d_dicom
    m = dcm_info.SharedFunctionalGroupsSequence.Item_1.PixelValueTransformationSequence.Item_1.RescaleSlope;
    b = dcm_info.SharedFunctionalGroupsSequence.Item_1.PixelValueTransformationSequence.Item_1.RescaleIntercept;
    % Note: the rescale intercept (b) is always 0 for PET images
else
    m = dcm_info.RescaleSlope;
    b = dcm_info.RescaleIntercept; % Note: the rescale intercept (b) is always 0 for PET images
end

data_in = m*data_in+b;

% get voxel dimensions
if is_3d_dicom
    pixel_size = dcm_info.SharedFunctionalGroupsSequence.Item_1.PixelMeasuresSequence.Item_1.PixelSpacing; % mm x mm
else
    pixel_size = dcm_info.(dicomlookup('0028','0030')); % mm x mm
end
vox_size = [pixel_size(1),pixel_size(2),dcm_info.(dicomlookup('0067','102F'))]; % mm x mm x mm

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
data_out = rot90(data_out,2);
fwrite(fopen(fname_amap,'w'),data_out,'int32');

%% write out phg_act_table

fid_phg_act_table = fopen(fname_phg_act_table,'w');
fprintf(fid_phg_act_table,'%d\n',n_bins);
fprintf(fid_phg_act_table,'%e\n',ac_res*(0:n_bins-1));
fclose('all');

end
