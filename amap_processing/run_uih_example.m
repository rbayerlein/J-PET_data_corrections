close all
clear
clc

%% user-configurable parameters

dir_in = 'cy3_test_101';
num_frames = 66;
frame_num = 61;

fname_amap = 'test_amap.raw';
fname_phg_act_table = 'test_act_table';

n_bins = 22100; % number of available bins

%% main code

data_out = SaveUih2Amap(dir_in,num_frames,frame_num,fname_amap,fname_phg_act_table,n_bins);
