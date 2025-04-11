function fname_out = GetFnameOut_umap(dir_in)
%GETFNAMEOUT_UMAP Summary of this function goes here
%   Detailed explanation goes here

if dir_in(end) == '/'
    dir_in = dir_in(1:end-1); % strip the / at the end for fileparts
end

[dir_out,basename_out] = fileparts(dir_in);
fname_out = strcat(dir_out,'/',basename_out,'.umap');

end
