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

%%%% STILL NEED TO ADD FILTER TO ELIMINATE MOVING CELLS - NEED PRACTICE SET
%%%%

%% INITIALIZE ALL PATHS

function main(tif_dir, out_dir, plate_name, file_name, fitc_str, dapi_str, size_cutoff, nuc_th, max_fitc, num_z, num_digit, num_img, P2A_min, max_col, max_row, max_sub, start_col, start_row, start_sub)
i_cnt = 0; % counter to track image number - should be 0, unless it crashed, then restart to # matching starting well
nuc_rad = sqrt(size_cutoff/pi);



%% READ ONE IMAGE FILE, CROP, ANALYZE AND SAVE EACH NUCLEUS (DAPI, MASK AND 3D FITC STACK)


% Reading image file from nikon microscope (i.e. one file for all channels, .nd2 extension).
% img  is a  1×4 cell array
% {3×2 cell}    {771 java.util.Hashtable}    {3×3 cell}    {1×1 loci.formats.ome.OMEPyramidStore}
% img{1} is num_ch×2 cell array
% e.g., 3 channels:
%     {2160×2560 uint16}    {'/data/mason_lab_imaging/DC_1gy_1hr_60x304.nd…'}
%     {2160×2560 uint16}    {'/data/mason_lab_imaging/DC_1gy_1hr_60x304.nd…'}
%     {2160×2560 uint16}    {'/data/mason_lab_imaging/DC_1gy_1hr_60x304.nd…'}

ch_nuc = 2; % channel for nucleus. Assume first one. Can be changed.
size_cutoff = 2000; % cutoff to remove small debris
nuc_th = 1.2;
nuc_rad = 40;
P2A_min = 2; % it's a high value compared to what I am used to, but for some reasons, many nuclei are deformed in these images
img = bfopen(image_path);
index_list=strfind(image_path,'/');
file_name = image_path(index_list(end)+1:end-4);
well_folder = image_path(1:index_list(end));
nuc_img = dip_image(img{1}{ch_nuc});
num_ch = size(img{1},1);
% segment image
mask = nuc_segmentor_local(nuc_img,nuc_rad,1,nuc_th,1,0);

% Measure number of nuclei
msr = measure(int16(mask),[],{'P2A','size'}); % label each cell
num_cells = length(msr); % number of cells is the length of measurement
final_mask = newim(size(mask));

% Save cropped mask and nuc and fitc
for i_nuc = 1:num_cells
    
    if (msr(i_nuc).P2A < P2A_min && msr(i_nuc).size>size_cutoff) % Only keep non deformed cells
        for i_ch = 1:num_ch
            [ch_img_cropped,crop_coordinates] = crop_from_mask(dip_image(img{1}{i_ch}),mask==i_nuc,10); % using 10 pixels to allow shift correction after acquisition , in case
            if i_ch == ch_nuc
                output_path = [well_folder, file_name,'DAPI_',int2str(i_nuc)];
                writeim(ch_img_cropped, output_path, 'TIFF',0); % actually writing file
                output_path = [well_folder, file_name,'MASK_',int2str(i_nuc)];
                writeim(mask(crop_coordinates(1,1):crop_coordinates(1,2),crop_coordinates(2,1):crop_coordinates(2,2)), output_path,'TIFF',0); % actually writing file
            end
            output_path = [well_folder, file_name,'_C',int2str(i_ch),'_',int2str(i_nuc)];
            writeim(ch_img_cropped, output_path,'ICSv1',0); % actually writing file
            final_mask = final_mask + (mask==i_nuc)*i_nuc;
        end
    end
end
% Save full image for visual inspection
mask = final_mask>0;
ofp = fopen([well_folder, file_name,'_nucleus_statistics.csv'],'w');
for i_ch = 1:num_ch
    msr = measure(int16(final_mask),dip_image(img{1}{i_ch}),{'Sum','Mean','StdDev','P2A','size'}); % measure basic statistics for nucleus in each channel
    if i_ch == 1
        msr_group = [msr.ID; msr.size; msr.P2A; msr.Mean; msr.Sum; msr.StdDev];
        msr_head = ["ID","Size","P2A",['Mean_Ch',int2str(i_ch)],['Sum_Ch',int2str(i_ch)],['StDev_Ch',int2str(i_ch)]];
    else
        msr_group = [msr_group; msr.P2A; msr.Mean; msr.Sum; msr.StdDev];
        msr_head = [msr_head, ['Mean_Ch',int2str(i_ch)],['Sum_Ch',int2str(i_ch)],['StDev_Ch',int2str(i_ch)]];
    end
end
textHeader = strjoin(msr_head,',');
fprintf(ofp,'%s\n',textHeader);
fclose(ofp);
dlmwrite([well_folder, file_name,'_nucleus_statistics.csv'],msr_group','append');

tmp_col = overlay(squeeze(stretch(clip(max(ch_img,[],3),max_fitc,200))),xor(mask,berosion(mask,2)),[255,0,0]);
well_info = sprintf('%s_%03d_color_FITC',well_name,sub_well); % creating rest of string for filename
output_path = fullfile(color_folder, well_info)
writeim(tmp_col, output_path, 'jpeg',0); % actually writing file
well_info = sprintf('%s_%03d_color_DAPI',well_name,sub_well); % creating rest of string for filename
tmp_col = overlay(squeeze(stretch(nuc_img)),xor(mask,berosion(mask,2)),[255,0,0]);
output_path = fullfile(color_folder, well_info)
writeim(tmp_col, output_path, 'jpeg',0); % actually writing file
fclose('all')
