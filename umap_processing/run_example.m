close all
clear
clc

%% user-configurable parameters

dir_in = {
    'cy1/CTAC_201';
    };

fname_out = cellfun(@GetFnameOut,dir_in,'UniformOutput',false);

hu_limits = [-1000,2000]; % lower and upper HU limits

fname_water = 'u_water_511keV.dat'; % 511 keV water density LUT

img_size_out = [239,239,679];
vox_size_out = 2.85*ones(1,3); % mm

tx_fov = 500; % mm

%% main code

figure();
for i = 1:numel(dir_in)
    [data_out,V] = Ct2Matmap(dir_in{i},hu_limits,fname_water,img_size_out,vox_size_out,tx_fov);
    fwrite(fopen(fname_out{i},'w'),data_out,'int32');
    subplot(numel(dir_in),1,i); histogram(V,'Normalization','probability'); title(dir_in{i},'Interpreter','none');
end
