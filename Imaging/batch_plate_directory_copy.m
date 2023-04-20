% script reads a sereies of plate directory and create another directory
% for each plate with only the jpg full view of each field and the
% processed jpg images with summary files. Useful to transfer data to
% collaborator w/o the huge size.
%
% S. Costes

% plate parameters

for i_plate = 313:314
    plate_name = sprintf('P%d',i_plate);
    target_dir_list = dir(['F:\Images\' plate_name '*']);
    num_plate_dir = length(target_dir_list);
    for i_dir = 1:num_plate_dir
        % Loop through all plate directory with same prefix P???
        if target_dir_list(i_dir).isdir
            target_dir =  ['F:\Images\' target_dir_list(i_dir).name '\'];
            copy_dir = ['F:\Images\' target_dir_list(i_dir).name '_small\'];
            fprintf('Directory %s contains data from %s\n',target_dir,plate_name);
            fprintf('Creating %s copy directory, smaller\n',copy_dir);
            mkdir(copy_dir);
            copyfile([target_dir '\*.txt'],copy_dir); % copy summary files
            file_list=dir(target_dir);
            ix = [file_list.isdir];
            file_list = file_list(ix);
            file_list = {file_list.name};
            num_file = length(file_list);
            num_list = 3:num_file;
            % start looping through all well directory in plate directory
            for i_list = num_list
                try
                    fprintf('Searching processed directory in %s\n',[target_dir file_list{i_list}]);
                    process_file = dir([target_dir file_list{i_list},'\Foci*']);
                    if ~isempty(process_file) % if more than one directory, only keep the good one
                        flag_good_dir = 1;
                        for i_pro = 1:length(process_file)
                            if flag_good_dir == 0 % if here, this means we have already a directory with good results. Get rid of the all the other ones
                                fprintf('Removing directory %s\n',[target_dir file_list{i_list} ,'\' process_file(i_pro).name]);
                                rmdir([target_dir file_list{i_list} ,'\' process_file(i_pro).name],'s');
                            else
                                process_dir = [target_dir file_list{i_list} ,'\' process_file(i_pro).name];
                                summary_file = dir([process_dir '\nuc_summary_foci_analysis*']);
                                if ~isempty(summary_file) % if there is a summary file, count number of nuclei
                                    fid = fopen([process_dir '\' summary_file.name]);
                                    datat = textscan(fid,'%f%s%f%f%f%f%f%f%s%f%f%f%f%*[^\n]','headerlines',1,'delimiter','\t');
                                    fclose(fid);
                                    num_nuc = length(datat{1});
                                    num_img = dir([process_dir '\*check.jpg']);
                                    num_img = length(num_img);
                                else % no file found, set the directory to be deleted
                                    num_nuc = 0;
                                    num_img = 1;
                                end
                                % Get rid of bad directories, if good stop
                                % searching other directories and set them
                                % to be deleted
                                if round(num_nuc/num_img) == 1
                                    flag_good_dir = 0;
                                    mkdir([copy_dir file_list{i_list}]); % create well directory
                                    copyfile([process_dir,'\*'],[copy_dir file_list{i_list}]); % copy processed jpg
                                    copyfile([target_dir file_list{i_list},'\*color_fitc*'],[copy_dir file_list{i_list}]); % copy fields
                                else
                                    fprintf('Removing directory %s\n',[target_dir file_list{i_list} ,'\' process_file(i_pro).name]);
                                    rmdir([target_dir file_list{i_list} ,'\' process_file(i_pro).name],'s');
                                end
                            end
                        end
                    end
                catch
                    fprintf('Error in directory %s\n',file_list{i_list});
                end
            end
        end
    end
end