close all
clear
clc

%% user-configurable parameters

dir_in = {
    '~/Desktop/human_test/CTAC_201'
    };

fname_out = cellfun(@GetFnameOut_umap,dir_in,'UniformOutput',false);

hu_limits = [-1000,2000]; % lower and upper HU limits

img_size_out = [239,239,679];
vox_size_out = 2.85*ones(1,3); % mm

tx_fov = 500; % mm

%% main code

for i = 1:numel(dir_in)
    [data_out,V] = Ct2Umap_only(dir_in{i},hu_limits,img_size_out,vox_size_out,tx_fov);
    fwrite(fopen(fname_out{i},'w'),data_out,'float');
end
