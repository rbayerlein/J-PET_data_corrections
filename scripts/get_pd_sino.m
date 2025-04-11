close all
clear
clc

%% user-configurable paramters

fname_prompt = 'block_sino_f0_prompts.raw';
fname_delay = 'block_sino_f0_randoms.raw';

fname_pd = 'block_sino_f0_pd.raw';

%% main code

data_prompt = fread(fopen(fname_prompt),'double');
data_delay = fread(fopen(fname_delay),'double');

data_pd = data_prompt - data_delay;

fwrite(fopen(fname_pd,'w'),data_pd,'double');

fclose('all');
