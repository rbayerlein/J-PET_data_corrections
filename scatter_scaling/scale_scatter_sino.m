% scatter scaling script for jPET data corrections
% using plane-by-plane scaling
% rbayerlein July 2024

% how to execute from terminal:
% matlab -nodesktop -nodisplay -r "scale_scatter_sino('${arg1}', '${arg2}', ...); quit"

function scale_scatter_sino(fname_prompt_sino, fname_rand_sino, fname_scatter_folder)

%% parameters
data_type_5D = 'uint16';
data_type_4D = 'double';
data_type_4D_r = 'single';

sino_size_5d = [51,39,25,25,20];
sino_size_4d = [51,39,25,25];

%% open files
disp("reading prompts sinogram");
p_sino5d = fread(fopen(fname_prompt_sino, 'rb'), inf, data_type_5D);
p_sino5d = reshape(p_sino5d, sino_size_5d);

disp("adding tof bins for prompts");
data_p = sum(p_sino5d,5);
        
disp("reading randoms sinogram");
r_sino4d = fread(fopen(fname_rand_sino, 'rb'), inf, data_type_4D_r);
r_sino4d = reshape(r_sino4d, sino_size_4d);

%%% TEMPORARY FIX FOR THE NUMBER OF RANDOMS IN The SINOGRAM %%%
r_sino4d = round(r_sino4d ./10);

disp("calculating prompts minus delayed sinogram");
pd_sino4d = data_p-r_sino4d;
pd_sino4d(pd_sino4d < 0 ) = 0; % catch negative entries

% free up space
clear r_sino4d p_sino5d data_p

disp("reading trues sinogram");
fname_t = [fname_scatter_folder, '/trues.sino5d'];
t_sino5d = fread(fopen(fname_t, 'rb'), inf, data_type_5D);
t_sino5d = reshape(t_sino5d, sino_size_5d);

disp("reading scatters sinogram");
fname_s = [fname_scatter_folder, '/scatters.sino5d'];
s_sino5d = fread(fopen(fname_s, 'rb'), inf, data_type_5D);
s_sino5d = reshape(s_sino5d, sino_size_5d);

disp("adding tof bins for trues");
data_t = sum(t_sino5d,5);

disp("adding tof bins for scatters");
data_s = sum(s_sino5d,5);

disp("calculating trues plus scatters");
ts_sino4d = data_t + data_s;

% free up space
clear data_t data_s t_sino5d

%% main function
disp("calculating scaling factors in plane-by-plane mode");
x_optimal = zeros(sino_size_4d(3),sino_size_4d(4));
for axA = 1 : sino_size_4d(3)
    for axB = 1 : sino_size_4d(4)
        sum_TS = sum(sum(ts_sino4d(:,:,axA,axB)));
        if sum_TS == 0 % catch division by zero
            x_optimal(axA,axB) = 0;
        else
            x_optimal(axA,axB) = sum(sum(pd_sino4d(:,:,axA,axB)))/sum_TS;
        end
        % apply non-zero negativity constraint to prevent errors
        if x_optimal(axA,axB) < 0
            x_optimal(axA,axB) = 0; 
        end
    end
end

disp("apply scaling to sinogram");
s_sino5d_scaled=zeros(sino_size_5d);
for axA = 1 : sino_size_4d(3)
    for axB = 1 : sino_size_4d(4)
        s_sino5d_scaled(:,:,axA,axB,:) = x_optimal(axA,axB) * s_sino5d(:,:,axA,axB,:);
    end
end


disp('applying smoothing to scaled sinos');
sm_kernel_FWHM = 4;    % in voxels
sm_sig = sm_kernel_FWHM/2.35;
for axB = 1:sino_size_4d(4)
    for axA = 1:sino_size_4d(3)
        for t = 1 : sino_size_5d(5)
            sum_before = sum(sum(s_sino5d_scaled(:,:,axA,axB,t)));
            if sum_before == 0
                continue
            end
            s_sino5d_scaled(:,:,axA,axB,t) = imgaussfilt(s_sino5d_scaled(:,:,axA,axB,t), sm_sig);
            sum_after = sum(sum(s_sino5d_scaled(:,:,axA,axB,t)));
            % correct for deviating number of counts
            s_sino5d_scaled(:,:,axA,axB,t) = s_sino5d_scaled(:,:,axA,axB,t)*(sum_before/sum_after);
        end
    end
end
clear s_sino5d

%% write out
disp("Saving data");
[FILEPATH,NAME,EXT] = fileparts(fname_s);
fname_s_scaled = strcat(FILEPATH, '/', NAME,'_scaled',EXT);
fname_scale_factor = strcat(FILEPATH,'/',NAME,'.scale_fac');
fprintf('Writing %s...\nand\n%s...\n',fname_s_scaled, fname_scale_factor);
fwrite(fopen(fname_s_scaled,'w'),s_sino5d_scaled, 'single');   % use single instead of double to limit file size.
fwrite(fopen(fname_scale_factor,'w'),x_optimal,'double');

end % function
