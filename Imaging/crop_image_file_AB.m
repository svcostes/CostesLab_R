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

tif_dir = ['/data/2022-11-15/895_Gamma_0_4h_Set3_A_21B'];
out_dir = ['/tmp/output/Gamma']; %Note used. Feel free to update path if you want things saved somewhere else
plate_name = ['895_Gamma_0_4h_Set3_A_21B'];
ch_name = {'C00(AF488)_M0000_ORG.tif','C03(AF594)_M0000_ORG.tif','C02(AF647)_M0000_ORG.tif'};
ch_name_short = {'AF488','AF594','AF647'};
ch_DAPI = 'C01(DAPI)';
Z_DAPI = 'Z0010'; % String determining which Z was saved for DAPI
num_ch = size(ch_name,2);
size_cutoff = 400; % cutoff to remove small debris
nuc_th = 4; % I had to set it high as the nuclei are super bright (400% above bdg)
nuc_rad = 15;
P2A_min = 1.5; % It will remove anything above 1.5 (i.e. non circular). You may have to increase it to keep less spherical nuclei
max_ch = [6000,6000,6000]; % Set maximum intensity for each channel in same order as ch_name.

%Reading image and setting file and folder names
image_path = [tif_dir,'/*',ch_DAPI,'*'];
D_img = dir(image_path); % List all files with DAPI string in it (assuming one DAPI per field)
num_field = size(D_img,1);

for i_field =1:num_field % Start reading each field and crop them
    image_path =[tif_dir,'/',D_img(i_field).name]; % first read nucleus field
    suffix_pos = findstr(Z_DAPI,image_path); % locate where the DAPI suffix starts
    field_name = image_path(1:suffix_pos-1);
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
            output_path = [field_name,'DAPI_',int2str(i_nuc)];
            [ch_img_cropped,crop_coordinates] = crop_from_mask(nuc_img,mask==i_nuc,5); % using 5 pixels to allow shift correction after acquisition , in case
            writeim(ch_img_cropped, output_path, 'TIFF',0); % actually writing file
            output_path = [field_name,'MASK_',int2str(i_nuc)];
            writeim(mask(crop_coordinates(1,1):crop_coordinates(1,2),crop_coordinates(2,1):crop_coordinates(2,2)), output_path,'TIFF',0); % actually writing file
            final_mask = final_mask + (mask==i_nuc)*nuc_cnt;
            nuc_cnt = nuc_cnt+1;
        end
    end
    num_cells = nuc_cnt -1;
    msr = measure(int16(final_mask),nuc_img,{'Sum','Mean','StdDev','P2A','size'}); %Get basic statistics for nucleus and DAPI intensity
    nuc_contour = xor(final_mask>0,bdilation(final_mask>0,3));
    tmp_col = overlay(squeeze(stretch(nuc_img)),nuc_contour,[255,0,0]);
    output_path = [field_name,'DAPI_2D'];
    writeim(tmp_col, output_path, 'jpeg',0); % actually writing file
    msr_group = [msr.ID; msr.size; msr.P2A; msr.Mean; msr.Sum; msr.StdDev];
    msr_head = ["ID","Size","P2A","Mean_DAPI","Sum_DAPI","StDev_DAPI"]; 
    
    %Second: read stacks for the other channels and save statistics
    for i_ch = 1:num_ch
        ch_img = readtimeseries([field_name,'Z*_',ch_name{i_ch}]); %Read channel image
        num_z = size(ch_img,3); % see how many Z
        mask = repmat(final_mask,1,1,num_z);
        msr = measure(int16(mask),ch_img,{'Sum','Mean','StdDev'}); % measure basic statistics for nucleus in each channel
        msr_group = [msr_group; msr.Mean; msr.Sum; msr.StdDev];
        msr_head = [msr_head, ['Mean_',ch_name_short{i_ch}],['Sum_',ch_name_short{i_ch}],['StDev_',ch_name_short{i_ch}]];
        for i_nuc = 1:num_cells
            [ch_img_cropped,crop_coordinates] = crop_from_mask(ch_img,mask==i_nuc,5); % using 5 pixels to allow shift correction after acquisition , in case
            output_path = [field_name,ch_name_short{i_ch},'_',int2str(i_nuc)];
            writeim(ch_img_cropped, output_path,'ICSv1',0); % actually writing file
        end
        tmp_col = overlay(squeeze(stretch(clip(max(ch_img,[],3),max_ch(i_ch),200))),nuc_contour,[255,0,0]);
        output_path = [field_name,ch_name_short{i_ch},'_2D'];
        writeim(tmp_col, output_path, 'jpeg',0); % actually writing file
    end
    
    % Save full image for visual inspection
    ofp = fopen([field_name,'nucleus_statistics.csv'],'w');
    textHeader = strjoin(msr_head,',');
    fprintf(ofp,'%s\n',textHeader);
    fclose(ofp);
    dlmwrite([field_name,'nucleus_statistics.csv'],msr_group','-append');
    fclose('all');
end
