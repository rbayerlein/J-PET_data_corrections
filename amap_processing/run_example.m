close all
clear
clc

%% user-configurable parameters

fname_in = 'lmrecon_explorer_OSEM_f60.intermediate.3';
fname_amap = 'test_amap.raw';
fname_phg_act_table = 'test_act_table';

n_bins = 22100; % number of available bins

img_size = [239,239,679];
vox_size = 2.85*ones(1,3); % mm

%% main code

data_out = SaveImg2Amap(fname_in,fname_amap,fname_phg_act_table,n_bins,img_size,vox_size);
