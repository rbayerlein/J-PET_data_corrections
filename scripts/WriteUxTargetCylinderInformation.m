function WriteUxTargetCylinderInformation(fid)
%WRITEUXTARGETCYLINDERINFORMATION Summary of this function goes here
%   Detailed explanation goes here

% *****THIS FUNCTION IS EXCLUSIVE TO THE J-PET SCANNER*****
% *****MODIFY AS NEEDED FOR OTHER SCANNERS*****

if nargin ~= 1
    error('Usage: WriteTargetUxCylinderInformation(fid)');
end

fprintf(fid,'%s\n','# TARGET CYLINDER INFORMATION');
fprintf(fid,'%s\n','NUM_ELEMENTS_IN_LIST target_cylinder = 3');
fprintf(fid,'\t%s\n','REAL target_zMin = -25.0');
fprintf(fid,'\t%s\n','REAL target_zMax = 25.0');
fprintf(fid,'\t%s\n','REAL radius = 33.936');
fprintf(fid,'%s\n','REAL acceptance_angle = 90.0');
fprintf(fid,'\n');

end
