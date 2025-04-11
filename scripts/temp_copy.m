amap_example.imgclose all
clear
clc

%% user inputs

user = 'rbayerlein';
server_name = 'exp-sim-001';

install_dir = '/home/rbayerlein/Code/Recon/scatter_correction/';

% dir_base = '/home/rbayerlein/Code/Recon/scatter_correction/dir_temp/test_phantom'; % base dir where intermediate files are saved
dir_base = '/mnt/data/rbayerlein/simset/test_phantom_reduced_nbins_1e10';

framing=[1,1320]; % num_frames, duration, [num_frames, duration, ...] - get from reconstruction GUI

% list mode file
% lm_file = '/media/rbayerlein/data/recon_data/lm_data_sets/toBeDeleted.lm'

% u-map
% dir_ctac = '/home/rbayerlein/Code/Recon/scatter_correction/recon_data/mumap/CTAC_201'; % get from reconstruction GUI
dir_ctac = [install_dir, 'dir_temp/recon_data/mumap/CTAC_201'];
tx_fov = 500; % mm

% a-map
num_iter = 2; % number of OSEM iterations - governs number of mcsc iterations
% dir_amaps_in = '/home/rbayerlein/Code/Recon/scatter_correction/recon_data/amap'; % must have the the same numel as number of frames
dir_amaps_in = [install_dir, 'dir_temp/recon_data/amap'];

% amap image dimensions (CT will be downscaled to match amap)
img_size = [239,239,679];
vox_size = [2.85,2.85,2.85]; % mm

% binaries paths
% bin_phg = '/home/rbayerlein/Code/Recon/scatter_correction/simset/2.9.2/bin/phg';
% bin_hist2lm = '/home/rbayerlein/Code/Recon/scatter_correction/simset/data_processing/hist2lm/bin/hist2lm';
% bin_lm2blocksino = '/home/rbayerlein/Code/Recon/scatter_correction/simset/data_processing/lm2blocksino_old/bin/lm2blocksino';
% bin_scatter_add_fac = '/home/rbayerlein/Code/Recon/scatter_correction/simset/data_processing/scatter_add_fac/bin/scatter_add_fac';

bin_phg = [install_dir, 'simset/2.9.2/bin/phg'];
bin_hist2lm = [install_dir, 'simset/data_processing/hist2lm/bin/hist2lm'];
bin_lm2blocksino =[install_dir, 'simset/data_processing/lm2blocksino_old/bin/lm2blocksino'];
bin_scatter_add_fac =[install_dir, 'simset/data_processing/scatter_add_fac/bin/scatter_add_fac'];

t12 = 109.771*60; % half-life (in s) - placeholder

% max number of threads to use per frame to match the reconstruction code
num_threads_per_frame = 8;
num_to_simulate = 1e10; % number of decays to simulate

% LUT
fname_lut = [bin_lm2blocksino(1:strfind(bin_lm2blocksino,'/bin')), 'include/index_blockpairs_transaxial_2x91x60_int16'];
%% setup

% create main folder structures
if not(isfolder(dir_base))
    CreateFolder(dir_base);
end
% create on sim node
cmd_0 = ['ssh ', user, '@', server_name, ' "mkdir ', dir_base,'"']; 
system(cmd_0);
pause(0.1);

dir_umap = strcat(dir_base,'/umap');
dir_amaps_out = strcat(dir_base,'/amap');
dir_phg_params = strcat(dir_base,'/phg_params');
dir_hist = strcat(dir_base,'/history');
dir_phg_log = strcat(dir_base,'/phg_log');

num_frames =1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

% cmd_rsync = ['rsync -arvu ', user, '@', server_name, ':', dir_hist, '/ ', dir_hist, '/ '];
% system(cmd_rsync);
% pause(0.1);
% cmd_rsync = ['rsync -arvu ', user, '@', server_name, ':', dir_umap, '/ ', dir_umap, '/ '];
% system(cmd_rsync);
% pause(0.1);
% cmd_rsync = ['rsync -arvu ', user, '@', server_name, ':', dir_amaps_out, '/ ', dir_amaps_out, '/ '];
% system(cmd_rsync);
% pause(0.1);
% cmd_rsync = ['rsync -arvu ', user, '@', server_name, ':', dir_phg_log, '/ ', dir_phg_log, '/ '];
% system(cmd_rsync);
% pause(0.1);
% cmd_rsync = ['rsync -arvu ', user, '@', server_name, ':', dir_phg_params, '/ ', dir_phg_params, '/ '];
% system(cmd_rsync);
% pause(0.1);

% convert history files to list-mode data

disp('Converting history files to .lm files...');
for i = 1:num_frames
    for j = 1:(num_iter-1) % first OSEM iteration is not scatter-corrected
        for k = 1:num_threads_per_frame
            cmd = sprintf('%s %s/f%05d/f%05d_part%02d.hist.%d %s/f%05d/f%05d_part%02d.%d',bin_hist2lm,dir_hist,i-1,i-1,k-1,j,dir_hist,i-1,i-1,k-1,j); % 0-index (i-1, k-1)
            system(cmd);
        end
    end
end

%% concatenate lm files

disp('Concatenating trues.lm files...');
for i = 1:num_frames
    for j = 1:(num_iter-1) % first OSEM iteration is not scatter-corrected
        cmd = sprintf('cat %s/f%05d/*.%d_trues.lm > %s/f%05d/f%05d.%d_trues.lm',dir_hist,i-1,j,dir_hist,i-1,i-1,j); % 0-index (i-1, k-1)
        system(cmd);
    end
end

disp('Concatenating scatters.lm files...');
for i = 1:num_frames
    for j = 1:(num_iter-1) % first OSEM iteration is not scatter-corrected
        cmd = sprintf('cat %s/f%05d/*.%d_scatters.lm > %s/f%05d/f%05d.%d_scatters.lm',dir_hist,i-1,j,dir_hist,i-1,i-1,j); % 0-index (i-1, k-1)
        system(cmd);
    end
end

%% convert lm to sino: run_lm2blocksino.sh

disp('Converting trues from list-mode to sinogram...');
for i = 1:num_frames
    for j = 1:(num_iter-1)
        cd(sprintf('%s/f%05d',dir_hist,i-1));
        cmd = sprintf('%s f%05d.%d_trues.lm f%05d.%d_trues.sino4d %s', bin_lm2blocksino, i-1, j, i-1, j, fname_lut)
        system(cmd);
    end
end

disp('Converting scatters from list-mode to sinogram...');
for i = 1:num_frames
    for j = 1:(num_iter-1)
        cd(sprintf('%s/f%05d',dir_hist,i-1));
        cmd = sprintf('%s f%05d.%d_scatters.lm f%05d.%d_scatters.sino4d %s', bin_lm2blocksino, i-1, j, i-1, j, fname_lut)
        system(cmd);
    end
end

system('echo $(date): All done!');