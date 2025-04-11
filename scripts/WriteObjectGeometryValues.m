function WriteObjectGeometryValues(fid,img_size,vox_size)
%WRITEOBJECTGEOMETRYVALUES Summary of this function goes here
%   Detailed explanation goes here

if nargin ~= 3
    error('Usage: WriteObjectGeometryValues(fid,img_size,vox_size)');
end

fprintf(fid,'%s\n','# OBJECT GEOMETRY VALUES');
fprintf(fid,'%s\n','BOOL point_source_voxels = false');
fprintf(fid,'%s\n','BOOL line_source_voxels = false');
fprintf(fid,'\n');
fprintf(fid,'%s%d\n','NUM_ELEMENTS_IN_LIST object = ',img_size(3)+1);
fprintf(fid,'\t%s%d\n','INT num_slices = ',img_size(3));

x_half = img_size(1)*vox_size(1)/2/10; % cm
y_half = img_size(2)*vox_size(2)/2/10; % cm
z_pos = -img_size(3)*vox_size(3)/2/10; % cm - iterates over slice number
for i0 = 0:img_size(3)-1 % 0-index
    fprintf(fid,'\t\t%s\n','NUM_ELEMENTS_IN_LIST slice = 9');
    fprintf(fid,'\t\t\t%s%d\n','INT slice_number = ',i0);
    fprintf(fid,'\t\t\t%s%f\n','REAL zMin = ',z_pos);
    z_pos = z_pos + vox_size(3)/10; % cm
    fprintf(fid,'\t\t\t%s%f\n','REAL zMax = ',z_pos);
    fprintf(fid,'\t\t\t%s%f\n','REAL yMin = ',-y_half);
    fprintf(fid,'\t\t\t%s%f\n','REAL yMax = ',y_half);
    fprintf(fid,'\t\t\t%s%f\n','REAL xMin = ',-x_half);
    fprintf(fid,'\t\t\t%s%f\n','REAL xMax = ',x_half);
    fprintf(fid,'\t\t\t%s%d\n','INT num_X_bins = ',img_size(1));
    fprintf(fid,'\t\t\t%s%d\n','INT num_Y_bins = ',img_size(2));
end
fprintf(fid,'\n');

end
