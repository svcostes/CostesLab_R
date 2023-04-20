% Create new directory structure mimicking high throughput output but
% reading data from CD7 (exported TIFF)
% This version creates a color_images folder, to allow QC after cropping.
% 
% Mason_crop_formt(tif_dir, out_dir, plate_name, file_name)
% tif_dir. Example: '/data/21B_TIFFS/2021-07-23'; % Where tiff from Zeiss CD7 have been put
%
% S. Costes, NASA, May 2019
% S. Costes, NASA, Sep 2019 - this version reads the direct tiff export
% from CD7 instead of the conversion from zvi to tiff.
% J. Oribello , 2021 - turn script into function on radbio
% Connie & Sylvain, 2021 - convert script to work for Mason DSUP data

%%%% STILL NEED TO ADD FILTER TO ELIMINATE MOVING CELLS - NEED PRACTICE SET
%%%%

%% INITIALIZE ALL PATHS

function main(tif_dir, out_dir, plate_name, file_name, fitc_str, dapi_str, size_cutoff, nuc_th, max_fitc, num_z, num_digit, num_img, P2A_min, max_col, max_row, max_sub, start_col, start_row, start_sub)
    i_cnt = 0; % counter to track image number - should be 0, unless it crashed, then restart to # matching starting well
    nuc_rad = sqrt(size_cutoff/pi);


    %%%%%%%!!!!!!!!!!!!!!!!!!!
    % This assumes that you always read a full export from one czi file. Which means plate was scanned on the first row from left to right
    col_increment = 1;
    %%%%%!!!!!!!!!!!!!!!!!!!!!

    %% READ IMAGE, CROP, ANALYZE AND SAVE
    crop_dir = fullfile(out_dir,plate_name); % fullfile for platform appropriate path creation
    fprintf('Saving Cropped images to "%s" \n',crop_dir);
    mkdir(crop_dir); % Create directory where cropped images are saved
    alphabet='A':'Z'; % to write well #
    log_file_path = fullfile(out_dir, plate_name, sprintf('%s_log.txt', plate_name));
    log_file = fopen(log_file_path, 'a');
    fprintf('Saving log to "%s" \n', log_file_path)

    % Reading is done by row first and then by column
    for i_col = start_col:max_col
        for i_row = start_row:max_row
            %start_row = 1; % once started, reseting to 1 to make sure when moving to next column, reading the full row, UNUSED in pipeline version
            well_name=[alphabet(i_row) num2str(i_col)]; % rename well to next one
            fprintf('Processing well %s\n',well_name);
            fprintf(log_file,'Processing well %s\n',well_name);
            well_folder = fullfile(out_dir, plate_name, well_name);
            color_folder = fullfile(out_dir, plate_name, 'Color_images');
            mkdir(well_folder);
            mkdir(color_folder);
            for sub_well = start_sub:max_sub
                % start_sub = 1; % once restarted, reseting to 1 to make sure when moving to next row, reading the full 15 sub-well, UNUSED in pipeline format
                fprintf('Well %s, Subwell#: %d\n',well_name,sub_well);
                fprintf(log_file,'Well %s, Subwell#: %d\n',well_name,sub_well);
                % read one image for DAPI
                dapi_name = fullfile(tif_dir, file_name, sprintf('%s_S*(P%d-%s)_Z%s_C01(%s)_M0000_ORG.tif',...
                    file_name,sub_well,well_name,int2str_format((num_z-1)/2,num_digit),dapi_str));
                dapi_read = dir(dapi_name);
                if length(dapi_read)~=1
                    fprintf('Problem reading file. %d image(s) for given well and position\nString search used:%s\n',...
                        length(dapi_read),dapi_name);
                    fprintf(log_file,'Problem reading file. %d image(s) for given well and position\nString search used:%s\n',...
                        length(dapi_read),dapi_name);
                    break;
                end
            image_path = fullfile(tif_dir, file_name, dapi_read.name)
                nuc_img = dip_image(imread(image_path));
                % read Z-stack for FITC
                fitc_img = newim(size(nuc_img,1),size(nuc_img,2),num_z);
                i_z = 0; i_zc = 0; % counter to add slices only if found. If not, print error, but continue to stacku up what you can
                while (i_zc<=num_z-1)
                    try
                        fitc_name = fullfile(tif_dir, file_name, sprintf('%s_S*(P%d-%s)_Z%s_C00(%s)_M0000_ORG.tif',...
                            file_name,sub_well,well_name,int2str_format(i_zc,num_digit),fitc_str));
                        fitc_read = dir(fitc_name);
                image_path = fullfile(tif_dir, file_name, fitc_read.name)
                        fitc_img(:,:,i_z) = dip_image(imread(image_path));
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
                    if (msr(i_nuc).P2A < P2A_min && msr(i_nuc).size>size_cutoff) % Only keep non deformed cells
                        [dapi_cropped,crop_coordinates] = crop_from_mask(nuc_img,mask==i_nuc,20); % using 10 pixels to allow shift correction after acquisition , in case
                        well_info = sprintf('%s_%03d_%03d_DAPI',well_name,sub_well,i_nuc); % creating rest of string for filename
                        output_path = fullfile(well_folder, well_info)
                        writeim(dapi_cropped, output_path, 'TIFF',0); % actually writing file
                        well_info = sprintf('%s_%03d_%03d_MASK',well_name,sub_well,i_nuc); % creating rest of string for filename
                        output_path = fullfile(well_folder, well_info)
                        writeim(mask(crop_coordinates(1,1):crop_coordinates(1,2),crop_coordinates(2,1):crop_coordinates(2,2)), output_path,'TIFF',0); % actually writing file
                        well_info = sprintf('%s_%03d_%03d_FITC',well_name,sub_well,i_nuc); % creating rest of string for filename
                        output_path = fullfile(well_folder, well_info)
                        writeim(fitc_img(crop_coordinates(1,1):crop_coordinates(1,2),crop_coordinates(2,1):crop_coordinates(2,2),:), output_path,'ICSv1',0); % actually writing file
                        final_mask = final_mask + (mask==i_nuc)*i_nuc;
                    end
                end
                % Save full image for visual inspection
                mask = final_mask>0;
                tmp_col = overlay(squeeze(stretch(clip(max(fitc_img,[],3),max_fitc,200))),xor(mask,berosion(mask,2)),[255,0,0]);
                well_info = sprintf('%s_%03d_color_FITC',well_name,sub_well); % creating rest of string for filename
                output_path = fullfile(color_folder, well_info)
                writeim(tmp_col, output_path, 'jpeg',0); % actually writing file
                well_info = sprintf('%s_%03d_color_DAPI',well_name,sub_well); % creating rest of string for filename
                tmp_col = overlay(squeeze(stretch(nuc_img)),xor(mask,berosion(mask,2)),[255,0,0]);
                output_path = fullfile(color_folder, well_info)
                writeim(tmp_col, output_path, 'jpeg',0); % actually writing file
            end
        end
    end
    fclose('all')
