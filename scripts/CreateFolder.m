function CreateFolder(dir_in)
%CREATEFOLDER Summary of this function goes here
%   Detailed explanation goes here

if exist(dir_in,'dir')
    error('%s exists!',dir_in);
else
    mkdir(dir_in);
end

end

