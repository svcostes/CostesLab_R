% Rename files and directory for cropped images due to error in creating
% them in the first place ;(
% It will only rename the files, not the directory. This should be done by
% hand to verify.
%
% S. Costes, NASA, May 2019

%% INITIALIZE ALL PATHS
global mic;
global shiftDAPI_FITC;
global microscope_portNumber;
global ASIDir;
global FIXED_INTENSITIES;
global max_num_nuc;
max_num_nuc = 200;
FIXED_INTENSITIES = 1;
shiftDAPI_FITC = -12;
microscope_portNumber = 5;
warning('OFF','MATLAB:nargchk:deprecated');
addpath('C:\Program Files\DIPimage 2.8\common\dipimage');
dip_initialise;
ASIDir = 'C:\Users\genuser\Documents\MATLAB\ASI\'; % Tiger direct PC
addpath(ASIDir);
addpath([ASIDir, '\work']);
addpath([ASIDir, '\work\script']);
addpath([ASIDir, '\work\Imaging']);
addpath([ASIDir, '\work\script\extern_communication']);
addpath([ASIDir, '\work\script\plates_loading']);
warning('off')


%% INPUT
out_dir = 'D:\'; % Directory where to save cropped images
plate_name = 'P457'; % This will name the new directory with cropped images and analysis.
num_digit = 4; % IF your numbers for position go until 99, use 2, if it goes until 999, use3. Default is 4 (1440 positions)
% Let us refer well in coordinate system. A1 is (1,1), H1 is (8,1)
first_well = [1,1];
last_well = [8,1];

%%%%%%%!!!!!!!!!!!!!!!!!!!
% This assumes that you always read a full export from one czi file. Which means plate was scanned on the first row from left to right
col_increment = 1;
%%%%%!!!!!!!!!!!!!!!!!!!!!

%% READ IMAGE, CROP, ANALYZE AND SAVE
crop_dir = ([tif_dir,'\',plate_name]);
mkdir(crop_dir); % Create directory where cropped images are saved
read_flag = 1; % set it to 0 once we get in last well
last_flag = 0; % flag it when reaching last well
i_cnt = 1; % counter to track image number
i_row = first_well(1,1);
i_col = first_well(1,2);

% find index for list of created folders sorted in time of creation
i_folder = 1; % increment tracking folder in order
list_dir = dir([out_dir '\' plate_name]);
valid_file = find([list_dir.isdir]);
valid_file = valid_file(3:end); % removing the . and .. directories
file_date = [list_dir.datenum];
[~,s_ind] = sort([file_date]);
s_ind = s_ind(ismember(s_ind,valid_file));
num_dir = length(s_ind(:));


max_col = 12;  % number of columns per plate
max_row = 8; % number of rows per plate
sub_well = 1; % variable used to tracks number of position per well. Once it reaches num_img, we move to next well
cell_count = 0; % variable counting # of cell in each well.
alphabet='A':'Z'; % to write well #
well_name=[alphabet(i_row) num2str(i_col)]; % initialize name of first well
fprintf('Processing well %s\n',well_name);

while read_flag
    if sub_well==1 && i_folder >25
        fprintf('Checking %s - should be %s\n',list_dir(s_ind(i_folder)).name,well_name);
        if ~strcmp(well_name,list_dir(s_ind(i_folder)).name) % need to rename directory and files
            cur_dir = [out_dir '\' plate_name '\' list_dir(s_ind(i_folder)).name];
            list_well = dir(cur_dir);
            valid_file = find(~[list_well.isdir]);
            for i_file = 1:length(valid_file(:))
                u_ind = findstr(list_well(valid_file(i_file)).name,'_');
                fprintf('movefile %s\\%s %s\\%s%s\n',cur_dir,list_well(valid_file(i_file)).name,cur_dir,well_name,list_well(valid_file(i_file)).name(u_ind(1):end));
                eval(sprintf('movefile %s\\%s %s\\%s%s',cur_dir,list_well(valid_file(i_file)).name,cur_dir,well_name,list_well(valid_file(i_file)).name(u_ind(1):end)));
            end
        end
    end
    i_cnt = i_cnt + 1;
    sub_well = sub_well + 1;
    if sub_well > num_img % we are done with this well, moving to next one. Renaming well.
        if last_flag % if we are already in last well and we are moving to next position. we are done. exit.
            read_flag = 0;
            break;
        end
        sub_well = mod(sub_well,num_img);
        i_col = i_col + col_increment;
        if (i_col > max_col) || (i_col <1) % We moved to next row. Renumber row. Keep col the same. Inverse col increment
            i_row = i_row + 1;
            col_increment = - col_increment;
            i_col = i_col + col_increment;
        end
        well_name=[alphabet(i_row) num2str(i_col)]; % rename well to next one
        i_folder = i_folder + 1;
        fprintf('Processing well %s\n',well_name);
        % check if reached last well. If so flag it
        if (i_row == last_well(1)) && (i_col == last_well(2))
            last_flag = 1;
        end
    end
    
end
fclose('all')