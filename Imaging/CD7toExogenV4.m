% Create new directory structure mimicking high throughput output but
% reading data from CD7 (exported TIFF)
% This version creates a color_images folder, to allow QC after cropping.
%
% S. Costes, NASA, May 2019

%%%% STILL NEED TO ADD FILTER TO ELIMINATE MOVING CELLS - NEED PRACTICE SET
%%%% 

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
addpath('C:\Program Files\DIPimage 2.9\common\dipimage');
dip_initialise;
ASIDir = 'C:\Users\Radbio Lab\Documents\MATLAB\ASI\'; % Tiger direct PC
addpath(ASIDir);
addpath([ASIDir, '\work']);
addpath([ASIDir, '\work\script']);
addpath([ASIDir, '\work\Imaging']);
addpath([ASIDir, '\work\script\extern_communication']);
addpath([ASIDir, '\work\script\plates_loading']);
warning('off')


%% INPUT
strain = 'human'; % set to 'mouse' or 'human' 
tif_dir = '\\Desktop-67r2rld\g-speed\BNL19A\TIFF-SetA'; % Where tiff from Zeiss CD7 have been put
fitc_str = 'AF488_ORG';  % Replace to the right suffix for Dapi images
dapi_str = 'DAPI_ORG';  % Replace with the right suffix for FITC images
out_dir = '\\Desktop-67r2rld\g-speed\BNL19A\Processed-SetA'; % Directory where to save cropped images
plate_name = 'P501'; % This will name the new directory with cropped images and analysis.
file_name = '501_Ar_3_4h_Set3A'
num_z = 21; % number of Z stack
num_digit = 4; % IF your numbers for position go until 99, use 2, if it goes until 999, use3. Default is 4 (1440 positions)
num_img = 15; % number of full images per well
% Let us refer well in coordinate system. A1 is (1,1), H1 is (8,1)
first_well = [1,1];
last_well = [7,12];
i_cnt = 1; % counter to track image number - should be 1, unless it crashed, then restart to # matching starting well

% set parameters specific to the imaging
switch strain
    case 'mouse' % skin mouse cells are larger. Need a larger cutoff
        size_cutoff = 1500; nuc_th = 1.1; max_fitc = 25000; % Maximum pixels intensity for fitc. This is used to diplay max projection only. It allows to compare different batches
    case 'human'
        size_cutoff= 450; nuc_th = 1.5; max_fitc = 15000;
end
nuc_rad = sqrt(size_cutoff/pi);


%%%%%%%!!!!!!!!!!!!!!!!!!!
% This assumes that you always read a full export from one czi file. Which means plate was scanned on the first row from left to right
col_increment = 1; 
%%%%%!!!!!!!!!!!!!!!!!!!!!

%% READ IMAGE, CROP, ANALYZE AND SAVE
crop_dir = ([out_dir,'\',plate_name]);
mkdir(crop_dir); % Create directory where cropped images are saved
read_flag = 1; % set it to 0 once we get in last well
last_flag = 0; % flag it when reaching last well
i_row = first_well(1,1);
i_col = first_well(1,2);
% if mod(i_row,2) > 0  % We start counting column positively on odd rows. Once reaching last column, we move to next row and increment negatively
%     col_increment = 1;
% else
%     col_increment = -1;
% end

max_col = 12;  % number of columns per plate
max_row = 8; % number of rows per plate
sub_well = 1; % variable used to tracks number of position per well. Once it reaches num_img, we move to next well
cell_count = 0; % variable counting # of cell in each well.
alphabet='A':'Z'; % to write well #
well_name=[alphabet(i_row) num2str(i_col)]; % initialize name of first well
well_folder = [out_dir '\' plate_name '\' well_name];
color_folder = [out_dir '\' plate_name  '\Color_images'];
mkdir(well_folder);
mkdir(color_folder);
log_file = fopen([out_dir '\' plate_name '\' plate_name 'log.txt'], 'a');
fprintf('Processing well %s\n',well_name);
fprintf(log_file,'Processing well %s\n',well_name);



while read_flag
    fprintf('Image #: %d, Subwell#: %d\n',i_cnt,sub_well);
    fprintf(log_file,'Image #: %d, Subwell#: %d\n',i_cnt,sub_well);
    % read one image for DAPI
    nuc_img = dip_image(imread(sprintf('%s\\%s\\%s_s%sz%s_%s.tif',tif_dir,file_name,file_name,int2str_format(i_cnt,num_digit),int2str_format((num_z-1)/2+1,2),dapi_str)));
    % read Z-stack for FITC
    fitc_img = newim(size(nuc_img,1),size(nuc_img,2),num_z);
    i_z = 1; i_zc = 1; % counter to add slices only if found. If not, print error, but continue to stacku up what you can
    while (i_zc<=num_z)
        try
        fitc_img(:,:,i_z-1) = dip_image(imread(sprintf('%s\\%s\\%s_s%sz%s_%s.tif',tif_dir,file_name,file_name,int2str_format(i_cnt,num_digit),int2str_format(i_zc,2),fitc_str)));
        i_z = i_z + 1;
        catch
            fprintf('%s\n Excluding slice:%s\n',lasterr,sprintf('%s\\%s\\%s_s%sz%s_%s.tif',tif_dir,file_name,file_name,int2str_format(i_cnt,num_digit),int2str_format(i_zc,2),fitc_str));
            fprintf(log_file,'%s\n Excluding slice:%s\n',lasterr,sprintf('%s\\%s\\%s_s%sz%s_%s.tif',tif_dir,file_name,file_name,int2str_format(i_cnt,num_digit),int2str_format(i_zc,2),fitc_str));
        end
        i_zc = i_zc + 1;
    end
    fitc_img = fitc_img(:,:,0:i_z-2);
    % segment image to identify nuclei
    mask = nuc_segmentor_local(nuc_img,nuc_rad,1,nuc_th,1,0);
            
    % Measure number of nuclei
    msr = measure(int16(mask),[],{'P2A','size'}); % label each cell
    num_cells = length(msr); % number of cells is the length of measurement
    final_mask = newim(size(mask));
    
    % Save cropped mask and nuc and fitc
    for i_nuc = 1:num_cells
        if (msr(i_nuc).P2A < 1.9 && msr(i_nuc).size>size_cutoff) % Only keep non deformed cells
            cell_count = cell_count + 1;
            [dapi_cropped,crop_coordinates] = crop_from_mask(nuc_img,mask==i_nuc,20); % using 10 pixels to allow shift correction after acquisition , in case
            well_info = sprintf('%s_%03d_%03d_DAPI',well_name,sub_well,i_nuc); % creating rest of string for filename
            writeim(dapi_cropped,[well_folder '\' well_info],'TIFF',0); % actually writing file
            well_info = sprintf('%s_%03d_%03d_MASK',well_name,sub_well,i_nuc); % creating rest of string for filename
            writeim(mask(crop_coordinates(1,1):crop_coordinates(1,2),crop_coordinates(2,1):crop_coordinates(2,2)),[well_folder '\' well_info],'TIFF',0); % actually writing file
            well_info = sprintf('%s_%03d_%03d_FITC',well_name,sub_well,i_nuc); % creating rest of string for filename
            writeim(fitc_img(crop_coordinates(1,1):crop_coordinates(1,2),crop_coordinates(2,1):crop_coordinates(2,2),:),[well_folder '\' well_info],'ICSv1',0); % actually writing file
            final_mask = final_mask + (mask==i_nuc)*i_nuc;
        end
    end
    % Save full image for visual inspection
    mask = final_mask>0;
    tmp_col = overlay(squeeze(stretch(clip(max(fitc_img,[],3),max_fitc,200))),xor(mask,berosion(mask,2)),[255,0,0]);
    well_info = sprintf('%s_%03d_color_FITC',well_name,sub_well); % creating rest of string for filename
    writeim(tmp_col,[color_folder '\' well_info],'jpeg',0); % actually writing file
    well_info = sprintf('%s_%03d_color_DAPI',well_name,sub_well); % creating rest of string for filename
    tmp_col = overlay(squeeze(stretch(nuc_img)),xor(mask,berosion(mask,2)),[255,0,0]);
    writeim(tmp_col,[color_folder '\' well_info],'jpeg',0); % actually writing file

    
    % increment counter and check if we moved to next well - Make sure if
    % we moved to next column all increments are fixed as well...
    i_cnt = i_cnt + 1;
    sub_well = sub_well + 1;
    if sub_well > num_img % we are done with this well, moving to next one. Renaming well.
        if last_flag % if we are already in last well and we are moving to next position. we are done. exit.
            read_flag = 0;
            fprintf('Number of cell counted for well %s: %d\n',well_name,cell_count);
            fprintf(log_file,'Number of cell counted for well %s: %d\n',well_name,cell_count);
            break;
        end
        sub_well = mod(sub_well,num_img);
        i_col = i_col + col_increment;
        if (i_col > max_col) || (i_col <1) % We moved to next row. Renumber row. Keep col the same. Inverse col increment
            i_row = i_row + 1;
            col_increment = - col_increment;
            i_col = i_col + col_increment;
        end
        fprintf('Number of cell counted for well %s: %d\n',well_name,cell_count);
        fprintf(log_file,'Number of cell counted for well %s: %d\n',well_name,cell_count);
        cell_count = 0;
        well_name=[alphabet(i_row) num2str(i_col)]; % rename well to next one
        well_folder = [out_dir '\' plate_name '\' well_name];
        mkdir(well_folder);
        fprintf('Processing well %s\n',well_name);
        fprintf(log_file,'Processing well %s\n',well_name);
        % check if reached last well. If so flag it
        if (i_row == last_well(1)) && (i_col == last_well(2))
            last_flag = 1;
        end
    end
    
end
fclose('all')