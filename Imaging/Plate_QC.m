% Read all color images full view and let user decide if need to keep or
% not.
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
addpath('C:\Program Files\DIPimage 2.9\common\dipimage');
dip_initialise;
warning('off')


%% INPUT
plate_dir = 'D:\BNL19A_cropped_plates\P609'; % Where tiff from Zeiss CD7 have been put
both_img = 0; % Flag to set to 1 if you want to see both DAPI and FITC for QC. Default FITC alone
%% Read images and let user decide
color_dir = [plate_dir '\Color_images\'];
% first check if file exist. If it does, resume where you left it.
try
    [i_list,well_list,sub_list,keep_list] = textread([plate_dir '\qc_file.txt'],'%d\t%s\t%s\t%d\n');
    i_init = 1441 %i_list(end)+1;
    qc_file = fopen([plate_dir '\qc_file.txt'], 'a');
    fprintf('Resuming QC to image %d\n',i_init);
catch
    i_init = 1; %read from start, since it was not ran before
    qc_file = fopen([plate_dir '\qc_file.txt'], 'w');
    fprintf('QC was never ran. Start from image 1\n');
end

list_file = dir([color_dir '*DAPI.jpg']);
[~,idx] = sort([list_file.datenum]);
num_file = length(list_file);

for i_img = i_init:num_file
    i_date = idx(i_img);
    if both_img>0
        dapi = readcolorim([color_dir, list_file(i_date).name]);
        fitc = readcolorim([color_dir, list_file(i_date).name(1:end-8), 'FITC.jpg']);
        comp_img = joinchannels('RGB',[fitc{1},dapi{1}],[fitc{2},dapi{2}],[fitc{3},dapi{3}]);
    else
        comp_img = readcolorim([color_dir, list_file(i_date).name(1:end-8), 'FITC.jpg']);
    end
    well_ind=findstr(list_file(i_date).name,'_');
    well_name = list_file(i_date).name(1:well_ind(1)-1);
    subwell_name = list_file(i_date).name(well_ind(1)+1:well_ind(2)-1);
    dipshow(10,comp_img);
    ButtonName = questdlg('Do you want to keep the image?', ...
        sprintf('Well %s - Image %s',well_name,subwell_name), ...
        'Yes', 'No', 'Quit','Yes');
    switch ButtonName
        case 'Yes'
            fprintf(qc_file,'%d\t%s\t%s\t1\n',i_img,well_name,subwell_name);
        case 'Quit'
            close(gcf);
            break;
        otherwise
            fprintf(qc_file,'%d\t%s\t%s\t0\n',i_img,well_name,subwell_name);
    end % switch
end
fclose('all')