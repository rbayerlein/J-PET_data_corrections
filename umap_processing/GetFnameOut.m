function fname_out = GetFnameOut(dir_in)
%GETFNAMEOUT Summary of this function goes here
%   Detailed explanation goes here

if dir_in(end) == '/'
    dir_in = dir_in(1:end-1); % strip the / at the end for fileparts
end

[dir_out,basename_out] = fileparts(dir_in);
fname_out = strcat(dir_out,'/',basename_out,'.matmap');

end
