% Create new directory structure mimicking high throughput output but
% reading data from CD7 (exported TIFF)
% This version creates a color_images folder, to allow QC after cropping.
%
% params_CD7toExogenV6(tif_dir, out_dir, plate_name, file_name)
% tif_dir. Example: '/data/21B_TIFFS/2021-07-23'; % Where tiff from Zeiss CD7 have been put
%
% S. Costes, NASA, May 2019
% S. Costes, NASA, Sep 2019 - this version reads the direct tiff export
% from CD7 instead of the conversion from zvi to tiff.
% J. Oribello , 2021 - turn script into function on radbio
% K. Rienecker 2022 - trying to work with the script to read four total
% S. Costes 2023 - removed the function and turned it back into a script
% channels



%% READ ONE IMAGE FILE, CROP, ANALYZE AND SAVE EACH NUCLEUS (DAPI, MASK AND 3D FITC STACK)


% Reading image file from nikon microscope (i.e. one file for all channels, .nd2 extension).
% img  is a  1×4 cell array
% {3×2 cell}    {771 java.util.Hashtable}    {3×3 cell}    {1×1 loci.formats.ome.OMEPyramidStore}
% img{1} is num_ch×2 cell array
% e.g., 3 channels:
%     {2160×2560 uint16}    {'/data/mason_lab_imaging/DC_1gy_1hr_60x304.nd…'}
%     {2160×2560 uint16}    {'/data/mason_lab_imaging/DC_1gy_1hr_60x304.nd…'}
%     {2160×2560 uint16}    {'/data/mason_lab_imaging/DC_1gy_1hr_60x304.nd…'}

parameter_file = '/home/genuser/local_repos/Imaging/parameter_files/counterstain_human_blood_params_CD7.txt';
params = load_params(parameter_file); % parameter for detection
tif_dir = params.tif_dir;
out_dir parameter_file = '/home/genuser/local_repos/Imaging/parameter_files/counterstain_human_blood_params_CD7.txt';
params = load_params(parameter_file); % parameter for detection
tif_dir = params.tif_dir;
out_dir = params.out_dir;
plate_name = params.plate_name;
ch_name = params.long_name; %name of each channel; e.g. C00(AF488)_M0000_ORG.tif
num_ch = params.num_ch;
ch_name_short = params.short_name; % short string for channel: e.g. AF488
Z_DAPI = params.Zdapi; % String determining which Z was saved for DAPI and identifies DAPI single slice
nuc_th = params.nuc_th;
nuc_rad = params.nuc_rad;
size_cutoff = pi*nuc_rad^2*0.75; % cutoff to remove small debris
P2A_min = params.p2a; % It will remove anything above 1.5 (i.e. non circular). You may have to increase it to keep less spherical nuclei
max_ch = params.max_stretch; % Set maximum intensity for each channel in same order as ch_name.
min_ch = params.min_stretch; % Set maximum intensity for each channel in same order as ch_name.
chan_str = params.chan_str; % Channel string used for the spot channel (%use for spot detection)

% Create directory to save outputs
mkdir(out_dir);
mkdir([out_dir,'/Color_images']);
mkdir([out_dir,'/Crop_images']);
%Reading image and setting file and folder names
image_path = [tif_dir,'/*',Z_DAPI,'*DAPI*'];
D_img = dir(image_path); % List all files with DAPI string in it (assuming one DAPI per field)
num_field = size(D_img,1);

% Create output stats structure header
msr_head = ["Field Name","Well","Position","ID","Size","P2A","Mean_DAPI","Sum_DAPI","StDev_DAPI"];
out_dir = params.out_dir;
plate_name = params.plate_name;
ch_name = params.long_name; %name of each channel; e.g. C00(AF488)_M0000_ORG.tif
num_ch = params.num_ch;
ch_name_short = params.short_name; % short string for channel: e.g. AF488
Z_DAPI = params.Zdapi; % String determining which Z was saved for DAPI and identifies DAPI single slice
nuc_th = params.nuc_th;
nuc_rad = params.nuc_rad;
size_cutoff = pi*nuc_rad^2*0.75; % cutoff to remove small debris
P2A_min = params.p2a; % It will remove anything above 1.5 (i.e. non circular). You may have to increase it to keep less spherical nuclei
max_ch = params.max_stretch; % Set maximum intensity for each channel in same order as ch_name.
min_ch = params.min_stretch; % Set maximum intensity for each channel in same order as ch_name.
chan_str = params.chan_str; % Channel string used for the spot channel (%use for spot detection)

% Create directory to save outputs
mkdir(out_dir);
mkdir([out_dir,'/Color_images']);
mkdir([out_dir,'/Crop_images']);
%Reading image and setting file and folder names
image_path = [tif_dir,'/*',Z_DAPI,'*DAPI*'];
D_img = dir(image_path); % List all files with DAPI string in it (assuming one DAPI per field)
num_field = size(D_img,1);

% Create output stats structure header
msr_head = ["Field Name","Well","Position","ID","Size","P2A","Mean_DAPI","Sum_DAPI","StDev_DAPI"];
for i_ch = 1:num_ch
    msr_head = [msr_head, ['Mean_',ch_name_short{i_ch}],['Sum_',ch_name_short{i_ch}],['StDev_',ch_name_short{i_ch}]];
end

for i_field =1:num_field % Start reading each field and crop them
    image_path =[tif_dir,'/',D_img(i_field).name]; % first read nucleus field
    suffix_pos = findstr(Z_DAPI,image_path); % locate where the DAPI suffix starts
    field_name = image_path(1:suffix_pos-1); % Field name with full path
    suffix_pos = findstr('/',field_name);
    field_name_short = field_name(suffix_pos(end)+1:end-1); % Field name for ouput file
    suffix_pos1 = findstr('(P',field_name_short);
    suffix_pos2 = findstr('-',field_name_short);
    position_number = str2num(field_name_short(suffix_pos1(end)+2:suffix_pos2(end)-1));
    well_name = field_name_short(suffix_pos2(end)+1:end-1);
    
    
    
    % Read nucleus image
    nuc_img = readim(image_path);
    
    %Since you saved all your images as tiff, I don't need bfopen
    % index_list=strfind(image_path,'/');
    % file_name = image_path(index_list(end)+1:end-4);
    % well_folder = image_path(1:index_list(end));
    % nuc_img = dip_image(img{1}{1});
    
    
    
    % segment image
    mask = nuc_segmentor_local(nuc_img,nuc_rad,1,nuc_th,1,0);
    
    % Measure number of nuclei
    msr = measure(int16(mask),[],{'P2A','size'}); % label each cell
    num_cells = length(msr); % number of cells is the length of measurement
    nuc_cnt = 1;
    final_mask = newim(size(mask));
    % First, remove too small nuclei and save cropped nucleus from DAPI 2D image
    for i_nuc = 1:num_cells
        if (msr(i_nuc).P2A < P2A_min && msr(i_nuc).size>size_cutoff) % Only keep non deformed cells
            output_path = [out_dir,'/Crop_images/',field_name_short,'_',int2str_format(nuc_cnt,2),'_DAPI'];
            [ch_img_cropped,crop_coordinates] = crop_from_mask(nuc_img,mask==i_nuc,5); % using 5 pixels to allow shift correction after acquisition , in case
            writeim(ch_img_cropped, output_path, 'TIFF',0); % actually writing file
            output_path = [out_dir,'/Crop_images/',field_name_short,'_',int2str_format(nuc_cnt,2),'_MASK'];
            writeim(mask(crop_coordinates(1,1):crop_coordinates(1,2),crop_coordinates(2,1):crop_coordinates(2,2)), output_path,'TIFF',0); % actually writing file
            final_mask = final_mask + (mask==i_nuc)*nuc_cnt;
            nuc_cnt = nuc_cnt+1;
        end
    end
    num_cells = nuc_cnt -1;
    if num_cells == 0
        fprintf('Field %s is empty\n',field_name_short);
    else
        msr = measure(int16(final_mask),nuc_img,{'Sum','Mean','StdDev','P2A','size'}); %Get basic statistics for nucleus and DAPI intensity
        nuc_contour = xor(final_mask>0,bdilation(final_mask>0,3));
        tmp_col = overlay(squeeze(stretch(nuc_img)),nuc_contour,[255,0,0]);
        output_path = [out_dir,'/Color_images/',field_name_short,'_DAPI_2D'];
        writeim(tmp_col, output_path, 'jpeg',0); % actually writing file
        msr_group = [i_field*ones(1,num_cells); msr.ID; msr.size; msr.P2A; msr.Mean; msr.Sum; msr.StdDev]; % output stats for nucleus. First field is the field number repeated for each nucleus
        
        %Second: read stacks for the other channels and save field statistics
        %inside msr)_group arrray
        for i_ch = 1:num_ch
            ch_img = readtimeseries([field_name,'Z*_',ch_name{i_ch}]); %Read channel image
            num_z = size(ch_img,3); % see how many Z
            mask = repmat(final_mask,1,1,num_z);
            msr = measure(int16(mask),ch_img,{'Sum','Mean','StdDev'}); % measure basic statistics for nucleus in each channel
            msr_group = [msr_group; msr.Mean; msr.Sum; msr.StdDev];
            for i_nuc = 1:num_cells
                [ch_img_cropped,crop_coordinates] = crop_from_mask(ch_img,mask==i_nuc,5); % using 5 pixels to allow shift correction after acquisition , in case
                output_path = [out_dir,'/Crop_images/',field_name_short,'_',int2str_format(i_nuc,2),'_',ch_name_short{i_ch}];
                writeim(ch_img_cropped, output_path,'ICSv1',0); % actually writing file
            end
            tmp_col = overlay(squeeze(stretch(clip(max(ch_img,[],3),max_ch(i_ch),min_ch(i_ch)))),nuc_contour,[255,0,0]);
            output_path = [out_dir,'/Color_images/',field_name_short,'_',ch_name_short{i_ch},'_2D'];
            writeim(tmp_col, output_path, 'jpeg',0); % actually writing file
        end
        
        % Add field statistics to overall statistics file and row headers
        if i_field == 1 % In case first field encountered, we need to define the larger group array that will have both nuc and channel stats (longer array)
            msr_group_cum = msr_group;
            field_header = string(repmat(field_name_short,num_cells,1));
            well_header = string(repmat(well_name,num_cells,1));
        else
            msr_group_cum = [msr_group_cum,msr_group];
            field_header = [field_header;string(repmat(field_name_short,num_cells,1))];
            well_header = [well_header;string(repmat(well_name,num_cells,1))];
        end
        % Save full image for visual inspection
    end
end

ofp = fopen([out_dir,'/Crop_images/',plate_name,'_nucleus_statistics.csv'],'w');
textHeader = strjoin(msr_head,',');
fprintf(ofp,'%s\n',textHeader);
num_row = size(field_header,1);
num_col = size(msr_group_cum,1);
for i_row=1:num_row
    fprintf(ofp,'%s,%s',field_header(i_row,:),well_header(i_row,:));
    for i_col = 1:num_col
        fprintf(ofp,',%f',msr_group_cum(i_col,i_row));
    end
    fprintf(ofp,'\n');
end
fclose(ofp);

% This is the part doing spot detection
analysis_foci_exogen_with_mask_wo_record_v2([out_dir,'/Crop_images/'], 1, params, chan_str, 1 )
