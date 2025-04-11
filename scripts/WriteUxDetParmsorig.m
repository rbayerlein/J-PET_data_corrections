function WriteUxDetParms(fid,fname_hist)
%WRITEUXDETPARMS Summary of this function goes here
%   Detailed explanation goes here

% *****THIS FUNCTION IS EXCLUSIVE TO THE uEXPLORER SCANNER*****
% *****MODIFY AS NEEDED FOR OTHER SCANNERS*****

if nargin ~= 2
    error('Usage: WriteDetParms(fid,fname_hist)');
end

fprintf(fid,'%s\n','ENUM detector_type = block');
fprintf(fid,'%s\n','BOOL do_forced_interaction = false');
fprintf(fid,'%s\n','REAL reference_energy_keV = 511');
fprintf(fid,'%s\n','REAL energy_resolution_percentage = 11.7');
fprintf(fid,'%s%s%s\n','STR history_file = "',fname_hist,'"');
fprintf(fid,'%s\n','STR history_params_file = "@simset/repo/history/scanner.histparms"');

n_units = 8; % number of PET units
% axial shift
ax_pos = -8.478750e+01; % cm - initialization
ax_shift = 24.225; % cm - initialization
fprintf(fid,'%s%d\n','NUM_ELEMENTS_IN_LIST blocktomo_num_rings = ',n_units);
for i0 = 0:n_units-1 % 0-index
    fprintf(fid,'%s%d\n','#Ring ',i0);
    fprintf(fid,'\t%s\n','NUM_ELEMENTS_IN_LIST blocktomo_ring_description_list = 3');
    fprintf(fid,'\t\t%s\n','STR blocktomo_ring_parameter_file = "@simset/repo/scanner/uexplorer/unit.ringparms"');
    fprintf(fid,'\t\t%s%f\n','REAL blocktomo_ring_axial_shift = ',ax_pos);
    fprintf(fid,'\t\t%s\n','REAL blocktomo_ring_transaxial_rotation = 0');
    ax_pos = ax_pos + ax_shift; % cm
end

end
