% Create new directory structure mimicking high throughput output but
% reading data from CD7 (exported TIFF)
% This version creates a color_images folder, to allow QC after cropping.
%
% S. Costes, NASA, May 2019
% S. Costes, NASA, Sep 2019 - this version reads the direct tiff export
% from CD7 instead of the conversion from zvi to tiff.
% This version does not crop. It only creates the full view in color. You
% can change max_fitc to see how the imaging intensity has changed. We
% initially set it to 15000, but as of today, 11/12/20, it seems much
% dimmer and we had to set it to 5000 to see clear fitc on the full field
% of view. This suggests, CD7 setting are different or staining is much
% weaker. 

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
tif_dir = 'H:\G-SPEED2\BNL_19C\19C TIFFs\Ar'; % Where tiff from Zeiss CD7 have been put
fitc_str = 'AF488';  % Replace to the right suffix for Dapi images
dapi_str = 'DAPI';  % Replace with the right suffix for FITC images
out_dir = '\\Radbio-dell\g-speed2\BNL_19C\19C Cropped and Processed\Test'; % Directory where to save cropped images
plate_name = 'P753'; % This will name the new directory with cropped images and analysis.
file_name = '753_Ar_0_4h_Set1_A_19C';
num_z = 21; % number of Z stack
num_digit = 4; % IF your numbers for position go until 99, use 2, if it goes until 999, use3. Default is 4 (1440 positions)
num_img = 15; % number of full images per well


i_cnt = 0; % counter to track image number - should be 0, unless it crashed, then restart to # matching starting well

% set parameters specific to the imaging
switch strain
    case 'mouse' % skin mouse cells are larger. Need a larger cutoff
        size_cutoff = 1500; nuc_th = 1.1; max_fitc = 25000; % Maximum pixels intensity for fitc. This is used to diplay max projection only. It allows to compare different batches
    case 'human'
        size_cutoff= 450; nuc_th = 1.5; max_fitc = 5000; % change this if not happy with color image. However, this should not happen unless imaging has changed
end
nuc_rad = sqrt(size_cutoff/pi);


%%%%%%%!!!!!!!!!!!!!!!!!!!
% This assumes that you always read a full export from one czi file. Which means plate was scanned on the first row from left to right
col_increment = 1;
%%%%%!!!!!!!!!!!!!!!!!!!!!

%% READ IMAGE, CROP, ANALYZE AND SAVE
crop_dir = ([out_dir,'\',plate_name]);
mkdir(crop_dir); % Create directory where cropped images are saved
max_col = 12;  % number of columns per plate
max_row = 8; % number of rows per plate
start_col = 1; % change this number to whatever column you want to start from
start_row = 1; % change this number to whatever row you want to start from
start_sub = 1; % change this number to restart from a specific subwell. % CHANGE HERE IF READING ONLY SPECIFIC SUBWELLS
alphabet='A':'Z'; % to write well #
log_file = fopen([out_dir '\' plate_name '\' plate_name 'log.txt'], 'a');

% Reading is done by row first and then by column
for i_col = start_col:max_col
    for i_row = start_row:max_row
        start_row = 1; % once started, reseting to 1 to make sure when moving to next column, reading the full row
        well_name=[alphabet(i_row) num2str(i_col)]; % rename well to next one
        fprintf('Processing well %s\n',well_name);
        fprintf(log_file,'Processing well %s\n',well_name);
        well_folder = [out_dir '\' plate_name '\' well_name];
        color_folder = [out_dir '\' plate_name  '\Color_images'];
        mkdir(well_folder);
        mkdir(color_folder);
        for sub_well = start_sub:15
            start_sub =1; % once restarted, reseting to 1 to make sure when moving to next row, reading the full 15 sub-well
             % CHANGE HERE IF READING ONLY SPECIFIC SUBWELLS
            fprintf('Well %s, Subwell#: %d\n',well_name,sub_well);
            fprintf(log_file,'Well %s, Subwell#: %d\n',well_name,sub_well);
            % read one image for DAPI
            dapi_name = sprintf('%s\\%s\\%s_S*(P%d-%s)_Z%s_C01(%s)_M0000_ORG.tif',...
                tif_dir,file_name,file_name,sub_well,well_name,int2str_format((num_z-1)/2,num_digit),dapi_str);
            dapi_read = dir(dapi_name);
            if length(dapi_read)~=1
                fprintf('Problem reading file. %d image(s) for given well and position\nString search used:%s\n',...
                    length(dapi_read),dapi_name);
                fprintf(log_file,'Problem reading file. %d image(s) for given well and position\nString search used:%s\n',...
                    length(dapi_read),dapi_name);
                break;
            end
            nuc_img = dip_image(imread([tif_dir,'\',file_name,'\',dapi_read.name]));
            % read Z-stack for FITC
            fitc_img = newim(size(nuc_img,1),size(nuc_img,2),num_z);
            i_z = 0; i_zc = 0; % counter to add slices only if found. If not, print error, but continue to stacku up what you can
            while (i_zc<=num_z-1)
                try
                    fitc_name = sprintf('%s\\%s\\%s_S*(P%d-%s)_Z%s_C00(%s)_M0000_ORG.tif',...
                        tif_dir,file_name,file_name,sub_well,well_name,int2str_format(i_zc,num_digit),fitc_str);
                    fitc_read = dir(fitc_name);
                    fitc_img(:,:,i_z) = dip_image(imread([tif_dir,'\',file_name,'\',fitc_read.name]));
                    i_z = i_z + 1;
                catch
                    fprintf('%s\n Excluding slice:%s\n',lasterr,fitc_read.name);
                    fprintf(log_file,'%s\n Excluding slice:%s\n',lasterr,fitc_read.name);
                end
                i_zc = i_zc + 1;
            end
            fitc_img = fitc_img(:,:,0:i_z-1);
            % segment image to identify nuclei
            mask = nuc_segmentor_local(nuc_img,nuc_rad,1,nuc_th,1,0);
            
            % Measure number of nuclei
            msr = measure(int16(mask),[],{'P2A','size'}); % label each cell
            num_cells = length(msr); % number of cells is the length of measurement
            final_mask = newim(size(mask));
            
            % Save cropped mask and nuc and fitc
            for i_nuc = 1:num_cells
                if (msr(i_nuc).P2A < 1.9 && msr(i_nuc).size>size_cutoff) % Only keep non deformed cells
%                     [dapi_cropped,crop_coordinates] = crop_from_mask(nuc_img,mask==i_nuc,20); % using 10 pixels to allow shift correction after acquisition , in case
%                     well_info = sprintf('%s_%03d_%03d_DAPI',well_name,sub_well,i_nuc); % creating rest of string for filename
%                     writeim(dapi_cropped,[well_folder '\' well_info],'TIFF',0); % actually writing file
%                     well_info = sprintf('%s_%03d_%03d_MASK',well_name,sub_well,i_nuc); % creating rest of string for filename
%                     writeim(mask(crop_coordinates(1,1):crop_coordinates(1,2),crop_coordinates(2,1):crop_coordinates(2,2)),[well_folder '\' well_info],'TIFF',0); % actually writing file
%                     well_info = sprintf('%s_%03d_%03d_FITC',well_name,sub_well,i_nuc); % creating rest of string for filename
%                     writeim(fitc_img(crop_coordinates(1,1):crop_coordinates(1,2),crop_coordinates(2,1):crop_coordinates(2,2),:),[well_folder '\' well_info],'ICSv1',0); % actually writing file
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
        end
    end
end
fclose('all')