function main(parameter_file, crop_dir, plate_num, Nmax, chan_str) 
%batch_plate_BNL18B
%  script reads a plate directory and process all images based on parameter
%  file. This function will make sure if directory has already been analyzed
%  that the outuput is not corrupted by matching the number of nuclei from
%  .data file with number of nuclei in summary file. If wrong, it will
%  delete the processed directory file and rerun the analysis for each
%  corrupted directory, then it will concatenate.
%
%ARGS:
%  parameter_file: str, path to parameter file used for spot detection
%  crop_dir: str, path to cropped plate directory. Contains well directories, A1,H7,etc.
%  plate_num: str, plate number, must match directory
%  Nmax: int, maximum number of nuclei to be processed for foci detection
%  chan_str: str, image channel to perform foci detection on

%check_args(["parameter_file","crop_dir","plate_num","Nmax","chan_str"])

params = load_params(parameter_file);
% filepath
% ROOT_LOCATION/PLATE_DIR/WELL_DIR(Also color images)
% crop_dir/plate_dir/well_dir
plate_name = sprintf('P%s',plate_num);
well_dir_list = dir(fullfile(crop_dir, '*'));
% remove linux "." and ".."
well_dir_list = well_dir_list(~ismember({well_dir_list.name},{'.','..','Color_images'}));
num_well_dir = length(well_dir_list);
for i_dir = 1:num_well_dir
  if well_dir_list(i_dir).isdir
    target_well_dir =  fullfile(crop_dir, well_dir_list(i_dir).name);
    fprintf('Directory %s contains data from %s\n',target_well_dir, plate_name);
    sprintf('Spot detecting: %s', target_well_dir)
    analysis_foci_exogen_with_mask_wo_record_v3(target_well_dir, 1, params, chan_str, Nmax, 1);
% 	  sprintf('Concatenating on %s', crop_dir)
  end
end
concatenate_well_Ames(plate_name,crop_dir);

end 
