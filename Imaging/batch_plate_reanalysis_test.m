% script reads a plate directory and process all images based on parameter
% file. This function will make sure if directory has already been analyzed
% that the outuput is not corrupted by matching the number of nuclei from
% .data file with number of nuclei in summary file. If wrong, it will
% deletec the processed directory file and rerun the analysis for each
% corrupted directory, then it will concatenate.

% plate parameters
parameter_file = 'C:\Users\Wolverine\Desktop\MATLAB\ASI\work\Imaging\human_blood_params.txt';
chan_str = 'FITC';
params = load_params(parameter_file);
dir_num = [362]; %[362 - updated on 9/30/18 by SM];
for i_plate = dir_num
    plate_name = sprintf('P%d',i_plate);
    target_dir_list = dir(['D:\' plate_name '*']);
    num_plate_dir = length(target_dir_list);
    for i_dir = 1:num_plate_dir
        % Loop throught all plate directory with same prefix P???
        if target_dir_list(i_dir).isdir
            target_dir =  ['D:\' target_dir_list(i_dir).name '\'];
            fprintf('Directory %s contains data from %s\n',target_dir,plate_name);
            file_list=dir(target_dir);
            ix = [file_list.isdir];
            file_list = file_list(ix);
            file_list = {file_list.name};
            num_file = length(file_list);
            num_list = 3:num_file;
            % start looping through all well directory in plate directory
            
          for i_list = num_list
              
                fprintf('Searching processed directory in %s\n',[target_dir file_list{i_list}]);
                        process_file = dir([target_dir file_list{i_list},'\Foci*']);
                        if ~isempty(process_file)
                            fprintf('Removing directory %s\n',[target_dir file_list{i_list} ,'\' process_file.name]);
                            rmdir([target_dir file_list{i_list} ,'\' process_file.name],'s');
                        end
                        analysis_foci_exogen_with_mask_wo_record([target_dir file_list{i_list}], 1, params, chan_str, 1);
        %               analysis_foci_exogen_laplace([target_dir file_list{i_list}], 1, params, 1);
               % catch
                fprintf('Error in directory %s\n',file_list{i_list});
          end
        end
    end
end

    

    