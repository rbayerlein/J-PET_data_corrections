function WriteFilePaths(fid,attenuation_indexes,activity_indexes,activity_table,detector_params_file)
%WRITEFILEPATHS Summary of this function goes here
%   Detailed explanation goes here

if nargin ~= 5
    error('Usage: WriteFilePaths(fid,attenuation_indexes,activity_indexes,activity_table,detector_params_file)');
end

fprintf(fid,'%s\n','# COHERENT ANGULAR DISTRIBUTION FILES');
fprintf(fid,'%s\n','STR coherent_scatter_table = "@simset/phg.data/phg_ad_files"');
fprintf(fid,'\n');
fprintf(fid,'%s\n','# ISOTOPE DATA');
fprintf(fid,'%s\n','STR isotope_data_file = "@simset/phg.data/isotope_positron_energy_data"');
fprintf(fid,'\n');
fprintf(fid,'%s\n','# ACTIVITY INDEX FILE');
fprintf(fid,'%s%s%s\n','STR activity_indexes = "',activity_indexes,'"');
fprintf(fid,'\n');
fprintf(fid,'%s\n','# ACTIVITY TABLE FILE');
fprintf(fid,'%s%s%s\n','STR activity_table = "',activity_table,'"');
fprintf(fid,'\n');
fprintf(fid,'%s\n','# ACTIVITY INDEX TO TABLE TRANSLATION FILE');
fprintf(fid,'%s\n','STR activity_index_trans = "@simset/phg.data/phg_act_index_trans"');
fprintf(fid,'\n');
fprintf(fid,'%s\n','# ATTENUATION INDEX FILE');
fprintf(fid,'%s%s%s\n','STR attenuation_indexes = "',attenuation_indexes,'"');
fprintf(fid,'\n');
fprintf(fid,'%s\n','# ATTENUATION TABLE FILE');
fprintf(fid,'%s\n','STR attenuation_table = "@simset/phg.data/phg_att_table"');
fprintf(fid,'\n');
fprintf(fid,'%s\n','# ATTENUATION INDEX TO TABLE TRANSLATION FILE');
fprintf(fid,'%s\n','STR attenuation_index_trans = "@simset/phg.data/phg_att_index_trans"');
fprintf(fid,'\n');
fprintf(fid,'%s\n','# DETECTOR PARAMETER FILE');
fprintf(fid,'%s%s%s\n','STR detector_params_file = "',detector_params_file,'"');

end
