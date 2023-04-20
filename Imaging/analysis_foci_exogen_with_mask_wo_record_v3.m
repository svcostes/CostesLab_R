function analysis_foci_exogen_with_mask_wo_record_v3(target_dir, single_nuc_flag, params, chan_str, Nmax, PRINT_FLAG )
% This function reads files in target_dir and analyse them using parameters
% stored inside test_params.txt file located in param_dir. All output files
% and image thumbnails are saved inside output_dir.
% It automatically segment nuclei and count foci for each of these nuclei,
% return mean intensity for each channel per nucleus as well to
% differentiate cell type.
% Usage: analysis_foci_function(target_dir,name_out,disp_size)
%
%USER INPUT INFO:
% target dir: directory path where images are located
% name_out: string used to save nuc_summary and foci_summary output files
% disp_size: % display of images. On linux, typically 85% (85)
%
% V2: version that only process field of view validated by QC.
% V3: version that stop processing once Nmax nuclei have been analyzed
% (SVC 10/16/2019)

if ~exist('PRINT_FLAG') || isempty(PRINT_FLAG); PRINT_FLAG = false; end
if ~exist('chan_str') || isempty(chan_str); chan_str = 'FITC'; end
warning('OFF','MATLAB:nargchk:deprecated'); % Get rid of warning due to old versoin of matlab medmad

output_dir = target_dir;

if ~exist('name_out')
    name_out = ['foci_analysis_' date];
end

if ~exist('disp_size')
    disp_size = 100; % Set it to 100 if you want full image thumbnail. 50 if you want 50% (recommended)
end

dipsetpref('NumberOfThreads',1)
if ~exist('target_dir')
    target_dir = uigetdir('','Please pick directory where images are stored');
end


% Get list of images
[file_list, params] = update_list(target_dir, params);

% Set some parameters constant
if ~isfield(params,'fitc_max_foci_size')
    params.fitc_max_foci_size = 1;
end
if ~isfield(params,'txred_max_foci_size')
    params.txred_max_foci_size = 1;
end
params.target_dir = target_dir;
params.output_dir = output_dir;
params.dxy = 6.5/40; % 40X pixel size in um
params.dz = 1;
pix_size = params.dxy^2;
params.txred_index = [];
dimen = params.dim; % if set to 2, then do everything on projections
flag_first = 1; % Flag used so that first image ready to have output print headers. Then set it to 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ANALYSIS:

% tag = num2str(round(datenum(datestr(clock))*100000));

%cd(output_dir)
foci_out_dir_name = ['Foci_Processed_', date];
out_dir = fullfile(output_dir, foci_out_dir_name); 

mkdir(out_dir);
log_file = fopen(fullfile(output_dir, sprintf('run_log_%s_%s.txt', date, name_out)), 'w');
param_file = fopen(fullfile(output_dir, sprintf('param_used_%s_%s.txt', date, name_out)), 'w');
struct_print(param_file,params,1);
fclose(param_file);
params

% Read QC file, if it exists - June 2019 - s. costes
plate_dir = fullfile(target_dir, '..');
% first check if file exist. If it does, resume where you left it.
try % change from textread to import to be able to read manual qc_file. Added another line to turn sub_list to 3 digits string 7/24/2019
%    [i_list,well_list,sub_list,keep_list] = textread([plate_dir '\Syl_qc_file.txt'],'%d\t%s\t%s\t%d\n');
    idata = importdata(fullfile(plate_dir, 'qc_file.txt'),'\t',0);
    well_list = idata.textdata(:,2);
    sub_list = idata.data(:,1);
    keep_list = idata.data(:,2);
    qc_flag = 1;
    fprintf('Found QC file %s - will use this file to decide which field is to be kept\n',[plate_dir '\qc_file.txt']);
catch
    fprintf('QC file was not found. Keeping everything.\n');
    qc_flag = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%READ AND PROCESS:
i = 1;
image_flag = 1;
while image_flag>0
    % Get updated list of images and modifiy num_file
    [file_list, params] = update_list(target_dir, params);
    num_file = length(file_list);
    % THIS IF STATEMENT SHOULD BE CHANGED TO ADD A CHECKED FILE IN
    % DIRECTORY TELLING SOFTWARE IF IMAGING IS COMPLETE
    if (i == num_file || i>Nmax) % if last image loaded, get out at end of loop
        image_flag=0;
    end
    try % This will write a log if image failed to process
        %% READING IMAGE
        % if QC flag exists, make sure the field read has been approved by
        % user - If cannot be found in QC file, by default, keep the field.
        if qc_flag>0
            well_ind=findstr(file_list{i},'_');
            well_name = file_list{i}(1:well_ind(1)-1);
            subwell_name = file_list{i}(well_ind(1)+1:well_ind(2)-1);
            % Replace string comparison with integer comparison for manual
            % QC file - 07/24/2019
            %sub_ind = find(and(strcmp(well_list,well_name),strcmp(sub_list,subwell_name)));
            sub_ind = find(and(strcmp(well_list,well_name),sub_list==eval(subwell_name)));
            if isempty(sub_ind)
                qc_pass = 1;
            else
                if keep_list(sub_ind)>0
                    qc_pass = 1;
                else
                    qc_pass = 0;
                end
            end
        else % just set parameter if no QC to always pass the QC test
            sub_ind=1;
            keep_list(1)=1;
            qc_pass = 1;
        end
        
        if ~(qc_pass) % If we can keep the field or the field was not QC, read images for the field
            fprintf('Well %s, SubWell %s was removed by QC\n',well_name,subwell_name);
        else
            nuc = readim(fullfile(target_dir, file_list{i}));
            try
                fitc_img = readim(fullfile(target_dir,[file_list{i}(1:end-8) chan_str '.ics']));
            catch
                % the file is not found - continue
                i = i+1;
                continue;
            end
            
            if dimen==2
                fitc_img = squeeze(max(fitc_img,[],3));
            end
            % This next line could be change by nuc_segmentor_local if we want
            % to optimize nuclear segmentation after acquisition. For now, it's
            % using the segmentation found during acquisition.
            mask = readim(fullfile(target_dir, [file_list{i}(1:end-8) 'MASK.tif']));
            % make sure mask is a label image with center being the mask
            if single_nuc_flag
                mask = label(mask == mask(fix(end/2),fix(end/2)));
                num_nuc = 1;
            else
                num_nuc = max(mask);
            end
            
            %%SEGMENT NUCLEI
            if PRINT_FLAG
                fprintf('dxy=%f,dz=%f\n',params.dxy,params.dz);
                fprintf('segmenting nuclei in %s\n',fullfile(target_dir, file_list{i}));
            end
            fprintf(log_file,'segmenting nuclei in %s\n',fullfile(target_dir, file_list{i}));
            mask = erosion(dip_image(int16(mask)),1);
            
            % added on new pc on 7/7/2014
            % REMOVED BY S. COSTES ON MAY 2019
            %         try % If the signal in middle of image too low and cells is not picked at all with fitc it will crash. Keep image in this case.
            %             mask_overlay = convhull(squeeze(and(mask,threshold(sum(fitc_img,[],3))))); % only keep overlaying mask. Wrap to remove holes
            %             if sum(mask_overlay>0)/sum(mask>0) < 0.5  % if large shift stop analysis for this cell
            %                 i = i+1;
            %                 continue;
            %             end
            %             mask = dip_image(int16(mask_overlay));
            %         end
            if dimen >2
                mask = repmat(mask,1,1,size(fitc_img,3));
                nuc = repmat(nuc,1,1,size(fitc_img,3));
            end
            
            %%PROCESS SPOT DETECTION
            
            %Initial nuc measurements
            try
                clear ms_nuc;
            end
            try
                clear ms_foci;
            end
            
            % Make sure if no txred, all parameters are still set to random val
            if isempty(params.txred_index)
                txred_img = newim(size(fitc_img));
                params.txred_max_stretch=40000;
                params.txred_min_stretch=1;
            end
            
            % Make sure if no fitc, all parameters are still set to random val
            if isempty(params.fitc_index)
                fitc_img = newim(size(txred_img));
                params.fitc_max_stretch=10000;
                params.fitc_min_stretch=1;
            end
            
            ms_nuc.nucID = 1:num_nuc;
            ms_nuc.image = repmat({file_list{i}(1:end-4)},1,num_nuc);
            ms_nuc.fitc_background = repmat(mean(fitc_img(mask==0)),1,num_nuc); % indicator of background outside nucleus
            ms_nuc.txred_background = repmat(mean(txred_img(mask==0)),1,num_nuc); % indicator of background outside nucleus
            if dimen>2
                ms_nuc_dapi = measure(squeeze(mask(:,:,0)),squeeze(nuc(:,:,0)),{'mean','size','p2a'}); % compute mean dapi for each nucleus
                ms_nuc_fitc = measure(squeeze(mask(:,:,0)),squeeze(mean(fitc_img,[],3)),{'mean','sum'});
                ms_nuc_txred = measure(squeeze(mask(:,:,0)),squeeze(mean(txred_img,[],3)),{'mean','sum'});
            else
                ms_nuc_dapi = measure(mask,nuc,{'mean','size','p2a'}); % compute mean dapi for each nucleus
                ms_nuc_fitc = measure(mask,fitc_img,{'mean','sum'});
                ms_nuc_txred = measure(mask,txred_img,{'mean','sum'});
            end
            ms_nuc.nuc_dapi = ms_nuc_dapi.mean;
            ms_nuc.nuc_area = ms_nuc_dapi.size*pix_size^2;
            ms_nuc.p2a = ms_nuc_dapi.p2a;
            ms_nuc.nfoci = zeros(1,num_nuc);
            ms_nuc.spot_mean = nan(1,num_nuc);
            ms_nuc.fitc_mean = ms_nuc_fitc.mean;
            ms_nuc.fitc_sum = ms_nuc_fitc.sum;
            ms_nuc.txred_mean = ms_nuc_txred.mean;
            ms_nuc.txred_sum = ms_nuc_txred.sum;
            ms_nuc.foci_dapi = nan(1,num_nuc);
            %Initialize spot_mask image. This is important in case no spot are
            %detected in the full image, as then the last spot image from the
            %previous iteration endup being reused giving artifacts
            spot_mask = newim(size(nuc));
            
            switch(params.spot_index) % Spot index indicates which channel to use as spot for to see their position relative to the other channel
                case 'fitc' % Main spots are in FITC
                    spot_img = fitc_img;
                    spot_min = params.fitc_min_stretch;
                    spot_max = params.fitc_max_stretch;
                otherwise
                    spot_img = txred_img;
                    spot_min = params.txred_min_stretch;
                    spot_max = params.txred_max_stretch;
            end
            if PRINT_FLAG
                fprintf('Spot detection\n');
            end
            [spot_struct,spot_mask,label_nuc,full_spot]= james_spot_detection7(spot_img, mask, params.max_foci_size, 100, params.min, params.k_val, 2, 1);
            if sum(((spot_mask>0)*mask))>0 % make sure at least one focus was detected somewhere to compute any distance.
                % label all foci within nucleus by the nuc #
                nuc_label_foci = dip_image(int16((spot_mask>0)*mask));
                ms_nuc_foci = measure(nuc_label_foci,spot_img,{'mean','size'}); % compute mean intensity of spot in its channel, and # of foci
                ms_foci_nucID = measure(spot_mask>0,mask,'mean'); % get nucID for each foci
                ms_foci_focidapi = measure(spot_mask>0,nuc,'mean'); % measure dapi intensity at foci location
                ms_foci_foci_int = measure(spot_mask>0,spot_img,'mean'); % measure foci intensity in spot channel
                [sorted_nucID,index_nucID] = sort(ms_foci_nucID.mean); % Get index for reordering foci by nuclear ID, not by position in image
                num_foci = size(ms_foci_nucID.ID,2);
                % Add image name and nucID to structure output for
                % reference and create final structure output
                ms_nuc.nfoci(ms_nuc_foci.ID) = ms_nuc_foci.size;
                ms_nuc.spot_mean(ms_nuc_foci.ID) = ms_nuc_foci.mean;
                ms_foci.nucID = ms_foci_nucID.mean(index_nucID);
                ms_foci.foci_dapi = ms_foci_focidapi.mean(index_nucID);
                ms_foci.spot_mean = ms_foci_foci_int.mean(index_nucID);
                ms_foci.image =repmat({file_list{i}(1:end-4)},1,num_foci);
                %Compute dapi mean at foci location
                ms_foci.nuc_dapi = [];
                for i_nuc = 1:num_nuc
                    ms_nuc.foci_dapi(i_nuc) = mean(ms_foci.foci_dapi(ms_foci.nucID==i_nuc));
                    ms_foci.nuc_dapi = [ms_foci.nuc_dapi repmat(ms_nuc.nuc_dapi(i_nuc),1,sum(ms_foci.nucID==i_nuc))];
                end
            else % Else for if statement for testing at least one focus was detected.
                num_foci = 0;
                ms_foci = [];
            end
            % WRITE OUTPUT FILE
            fprintf(log_file,'Found %d nuclei and %d foci.\n',num_nuc,num_foci);
            if (flag_first == 1) % print header with output if first time run
                nuc_sum_fh = fopen(fullfile(out_dir,  ['nuc_summary_' name_out '.txt']), 'w');
                foci_sum_fh = fopen(fullfile(out_dir, ['foci_summary_' name_out '.txt']), 'w');
                struct_print(nuc_sum_fh,ms_nuc,1);
                struct_print(foci_sum_fh,ms_foci,1);
                flag_first = 0;
            else
                nuc_sum_fh = fopen(fullfile(out_dir, ['nuc_summary_' name_out '.txt']), 'at');
                foci_sum_fh = fopen(fullfile(out_dir, ['foci_summary_' name_out '.txt']), 'at');
                struct_print(nuc_sum_fh,ms_nuc,0);
                struct_print(foci_sum_fh,ms_foci,0);
            end
            fclose(nuc_sum_fh);
            fclose(foci_sum_fh);
            
            % WRITE OUTPUT IMAGES
            %
            % Prepare output image
            % Project 3D images on 2D plane
            if size(nuc,3) > 1
                spot2D = squeeze(max(spot_img,[],3));
                mask2D = squeeze(max(mask,[],3));
                nuc2D = squeeze(max(nuc,[],3));
                spot_mask2D = squeeze(max(spot_mask,[],3));
            else
                spot2D = spot_img;
                mask2D = mask;
                nuc2D = nuc;
                spot_mask2D = spot_mask;
            end
            
            
            %CREATE SPOT 2D DISPLAY
            spot2D = stretch(spot2D, 0, 100, min(spot2D)/spot_min*10,max(spot2D)/spot_max*255); % Stretch all fitc images the same way
            spot2D(spot2D>255) = 255; % this is to avoid weird cliping of image above 255 intensity with writeim function
            overlay_img = overlay_phase(spot2D,stretch(nuc2D),'blue');
            overlay_img = overlay(overlay_img,xor(bdilation(spot_mask2D>0,2),bdilation(spot_mask2D>0)),[255 0 0]);
            edge = norm(gradient(mask2D))>0;
            edge = bskeleton(edge);
            overlay_img = overlay(overlay_img,edge,[255 255 0]);
            writeim(overlay_img,fullfile(out_dir, [file_list{i}(1:end-4) '_spot_check.jpg']),'jpeg');
            overlay_img = overlay(spot2D,edge,[255 0 0]);
            writeim(overlay_img,fullfile(out_dir, [file_list{i}(1:end-4) '_spot_check_no_foci_contour.jpg']),'jpeg');
        end
    catch MError
        stack=MError.stack;
        file_err=stack.file;
        name_err=stack.name;
        line_err=stack.line;
        error_message=[MError.identifier 10 MError.message 10 file_err 10 name_err 10 'line ' num2str(line_err)];
        fprintf(log_file, 'Could not execute analysis_foci_exogen_with_mask\n%s\n',error_message);
        rethrow(MError);
    end
    i = i + 1; % increment file count
end

fclose(log_file);

function [file_list, params] = update_list(target_dir, params)

%get rid of directories in listing
file_list = dir(fullfile(target_dir, '*DAPI.tif')); %DAPI listed as tif on ASI
if isempty(file_list)
    throw(MException('Id:EmptyDir', 'No DAPI images in this directory'));
end
[date_s,date_ind] = sort(cell2mat({file_list.datenum}));
file_list = {file_list.name};
file_list = file_list(date_ind);

%for output
dot_index = findstr(file_list{1}, '.');
params.suffix = file_list{1}(dot_index+1:end);
params.range_low = 1;
params.range_high = length(file_list);
