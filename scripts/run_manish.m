close all
clear
clc

%% user inputs
% define base directory where temporary output should be saved
dir_base = '/data/4/users/mdas/Softwares/simset/iq_f18_test_new';

% u-map
fname_umap = '/data/4/users/mdas/iq_f18/simset_input/test_umap_resampled.img';
tx_fov = 509.147; % mm

% activity map

dir_amaps_in = '/data/4/users/mdas/iq_f18/simset_input/';

% amap and umap image dimensions
img_size = [200,200,200];
vox_size = [2.5,2.5,2.5]; % mm

% binaries paths
simset_path = '/data/4/users/mdas/Softwares/simset/';

bin_phg = [simset_path,'/2.9.2_new_mat_table/bin/phg'];
bin_hist2lm = [simset_path, '/data_processing/hist2lm/bin/hist2lm'];
bin_lm2blocksino = [simset_path,'/data_processing/lm2blocksino_old/bin/lm2blocksino'];
bin_scatter_add_fac = [simset_path, '/data_processing/scatter_add_fac/bin/scatter_add_fac'];

% max number of threads to use per frame to match the reconstruction code
num_threads_per_frame = 30;
num_to_simulate = 1000000000 ; %5e9;%10e10; % number of decays to simulate in total
% the decays will be equally distributed over all threads.

%% setup

% create main folder structures
if not(isfolder(dir_base))
    CreateFolder(dir_base);
end

dir_umap = strcat(dir_base,'/umap');
dir_amaps_out = strcat(dir_base,'/amap');
dir_phg_params = strcat(dir_base,'/phg_params');
dir_hist = strcat(dir_base,'/history');
dir_phg_log = strcat(dir_base,'/phg_log');

CreateFolder(dir_umap);
CreateFolder(dir_amaps_out);
CreateFolder(dir_phg_params);
CreateFolder(dir_hist);
CreateFolder(dir_phg_log);

% u-map
hu_limits = [-1000,2000]; % lower and upper HU limits
fname_water = 'u_water_511keV.dat'; % 511 keV water density LUT

% a-map
n_bins = 2301; % number of available bins - same as number of materials


%% starting

system('echo $(date): Starting run_example.m...');

%% generate and save material index map


% % start parallel pool
parpool(maxNumCompThreads-1); % leave 2 threads out to minimize SSH lag
 
disp('Generating material map...');
 
cd ../umap_processing/

fname_matmap = strcat(dir_umap,'/CTAC_201.matmap');
data_matmap = Mumap2Matmap_JPET(fname_umap, fname_water,img_size);
fwrite(fopen(fname_matmap,'w'),data_matmap,'int32');
cd ../scripts/

% shut down parallel pool
delete(gcp('nocreate'))

%% generate and save each activity map (w/ corresponding activity table) and phg_params

disp('Generating activity map(s) and phg_params...');

rng('shuffle'); % initialize RNG
num_to_simulate_per_thread = num_to_simulate/num_threads_per_frame;

cd ../amap_processing/
fname_amaps_in = [dir_amaps_in, 'recon_norm_rand_it4.img'];
fname_amaps_out = sprintf('%s/amap',dir_amaps_out);
fname_phg_act_tables = sprintf('%s/phg_act_table',dir_amaps_out); 
SaveImg2Amap_original(fname_amaps_in,fname_amaps_out,fname_phg_act_tables,n_bins,img_size,vox_size);
cd ../scripts/

% create a phg_params file for each chopped frame
for k = 1:num_threads_per_frame

    % generate and save each *.detparms

    fname_detector_params = sprintf('%s/part%02d.detparms',dir_phg_params, k-1);
    fid_detector_params = fopen(fname_detector_params,'w');

    fname_hist = sprintf('%s/part%02d.hist',dir_hist,k-1);
    WriteUxDetParms(fid_detector_params,fname_hist);
    fclose(fid_detector_params);

    % generate and save each *.phg_params

    fname_phg_params = sprintf('%s/part%02d.phg_params',dir_phg_params,k-1); 
    fid_phg_params = fopen(fname_phg_params,'w');
    x_rng = randi(1e6); % random seed
    WritePhgOptions(fid_phg_params,num_to_simulate_per_thread,x_rng);
    WriteObjectGeometryValues(fid_phg_params,img_size,vox_size);
    WriteUxTargetCylinderInformation(fid_phg_params);
    WriteFilePaths(fid_phg_params,fname_matmap,fname_amaps_out,fname_phg_act_tables,fname_detector_params);
    fclose(fid_phg_params);

end

%% setting up simulation scripts

% create Bash scripts for each frame
disp('Setting up simulation scripts...');

fname_phg_sh = sprintf('%s/run.sh',dir_phg_params);
fid_phg_sh = fopen(fname_phg_sh,'w');
fprintf(fid_phg_sh,'%s\n','#!/bin/bash');
fprintf(fid_phg_sh,'\n');
fprintf(fid_phg_sh,'%s%d\n','num_threads=',num_threads_per_frame);
fprintf(fid_phg_sh,'\n');
fprintf(fid_phg_sh,'%s%s\n','fname_simset=',bin_phg);
fprintf(fid_phg_sh,'\n');
fprintf(fid_phg_sh,'%s%s\n','dir_phg_params=',dir_phg_params);
fprintf(fid_phg_sh,'%s%s\n','dir_simset_log=',dir_phg_log);
fprintf(fid_phg_sh,'\n');
fprintf(fid_phg_sh,'%s\n','echo "$(date): Starting simulations..."');
fprintf(fid_phg_sh,'\n');
fprintf(fid_phg_sh,'%s\n','parallel --will-cite -j ${num_threads} "${fname_simset} ${dir_phg_params}/{} > ${dir_simset_log}/{}.log" ::: *.phg_params');
fprintf(fid_phg_sh,'\n');
fprintf(fid_phg_sh,'%s\n','echo "$(date): Finished simulations."');
fclose(fid_phg_sh);



%% run simulations

disp('Starting simulations. This may take a while...');
cd(sprintf('%s',dir_phg_params)); % script must be run in working directory
fname_phg_sh = sprintf('%s/run.sh',dir_phg_params);
cmd = sprintf('%s%s','source ',fname_phg_sh);
system(cmd);
    

%% convert history files to list-mode data

%disp('Converting history files to .lm files...');

%for k = 1:num_threads_per_frame
    %cmd = sprintf('%s %s/part%02d.hist %s/part%02d',bin_hist2lm,dir_hist,k-1,dir_hist,k-1); 
    %system(cmd);
%end


disp('Converting history files to .lm files...');

% Prepare the command template for GNU Parallel
cmd_template = '%s %s/part%02d.hist %s/part%02d';
cmd_base = bin_hist2lm; 

% Prepare commands for GNU Parallel
commands = cell(1, num_threads_per_frame);
for k = 1:num_threads_per_frame
    cmd = sprintf(cmd_template, cmd_base, dir_hist, k-1, dir_hist, k-1);
    commands{k} = cmd;
end

% Convert MATLAB cell array to a single string separated by newlines
cmd_str = strjoin(commands, '\n');

% Use GNU Parallel to execute the commands in parallel
gnu_parallel_cmd = ['echo -e "', cmd_str, '" | parallel --eta'];

% Execute GNU Parallel command using system()
system(gnu_parallel_cmd);

disp('Conversion completed.');




%% concatenate lm files

disp('Concatenating trues.lm files...');
cmd = sprintf('cat %s/*_trues.lm > %s/trues.lm',dir_hist,dir_hist);
system(cmd);

disp('Concatenating scatters.lm files...');
cmd = sprintf('cat %s/*_scatters.lm > %s/scatters.lm',dir_hist,dir_hist);
system(cmd);

fInfo = dir(sprintf('%s/trues.lm', dir_hist));
fileSize_t = fInfo(1).bytes;
fInfo = dir(sprintf('%s/scatters.lm', dir_hist));
fileSize_s = fInfo(1).bytes;

SF = fileSize_s/(fileSize_s+fileSize_t)*100;
fprintf('Scatter fraction: %2.2f %%\n', SF);

%% all done!
system('echo $(date): All done!');
