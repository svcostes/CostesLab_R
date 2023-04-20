% script reads a plate directory and process all images based on parameter
% file

target_dir = 'F:\Images\P303_73676278850\';
parameter_file = 'C:\Users\genuser\Documents\MATLAB\Microscope-Control\work\Imaging\mouse_skin_params.txt';
% target_dir = 'F:\Images\P0095_73635180838\';
% parameter_file = 'C:\Users\genuser\Documents\MATLAB\Microscope-Control\work\Imaging\mouse_params.txt';
file_list=dir(target_dir);
ix = [file_list.isdir];
file_list = file_list(ix);
file_list = {file_list.name};
num_file = length(file_list);
params = load_params(parameter_file);
chan_str = 'FITC';
num_list = 3:num_file;
num_list = [12,13,14,22,24,32,34]-1;


for i_list = num_list

    fprintf('Searching processed directory in %s\n',[target_dir file_list{i_list}]);
        process_file = dir([target_dir file_list{i_list},'\Foci*']);
        if ~isempty(process_file)
            fprintf('Removing directory %s\n',[target_dir file_list{i_list} ,'\' process_file.name]);
            rmdir([target_dir file_list{i_list} ,'\' process_file.name],'s');
        end
        analysis_foci_exogen_with_mask_wo_record([target_dir file_list{i_list}], 1, params, chan_str, 1);
        %   analysis_foci_exogen_laplace([target_dir file_list{i_list}], 1, params, 1);
      catch
        fprintf('Error in directory %s\n',file_list{i_list});
    end
end
