function [record] = analysis_foci_exogen_with_mask_v2(target_dir, record, well_struct, single_nuc_flag, PRINT_FLAG )
% This function reads files in target_dir and analyse them using parameters
% based on stored parameters for skin or blood, mice or humans, 53bp1 or h2ax or both.
% All output files
% and image thumbnails are saved inside output_dir.
% It automatically segment nuclei and count foci for each of these nuclei,
% return mean intensity for each channel per nucleus as well to
% differentiate cell type.
% Usage: analysis_foci_exogen_with_mask_v2(target_dir, record, well_struct, single_nuc_flag, PRINT_FLAG )
%
%USER INPUT INFO:
% target dir: directory path where images are located
% name_out: string used to save nuc_summary and foci_summary output files
% disp_size: % display of images. On linux, typically 85% (85)
%
% This version should be able to handle two colors at once.
% S.V. Costes, September 2016

if ~exist('PRINT_FLAG') || isempty(PRINT_FLAG); PRINT_FLAG = false; end

print_label = 0; % make it positive if you want to save images with label #
% Printing label requires to display each image on the
% screen (slow)
output_dir = target_dir;
disp_size = 100; % Set it to 100 if you want full image thumbnail. 50 if you want 50% (recommended)

if isunix % Determine slash direction
    slash_str = '/';
else
    slash_str = '\';
end

dipsetpref('NumberOfThreads',1)
if ~exist('target_dir')
    target_dir = uigetdir('','Please pick directory where images are stored');
end

% Make sure path are right for unix
if isunix %if unix, make sure all slash are in the right direction or it wont find the file...
    rep_ind = findstr(target_dir,'\');
    for i=1:length(rep_ind)
        target_dir(rep_ind) = '/';
    end
end

% Make sure only one slash at the end of directory path
if ~strcmp(slash_str,target_dir(end))
    target_dir = [target_dir,slash_str];
end

% loading params
switch lower(well_struct.Race)
    case 'mouse'
        switch lower(well_struct.Cell_0x20_Type)
            case 'PBMC'
                fprintf('Loading mouse blood parameters\n');
                params = load_params('all_mouse_blood_params.txt');
            otherwise
                fprintf('Loading mouse skin parameters for 53BP1\n');
                params = load_params('all_mouse_skin_params.txt');
        end
    otherwise % default human
        fprintf('Loading human blood parameters ');
        params = load_params('all_human_blood_params.txt');
end

% show parameters...
fprintf('Processing %s, with following parameters settings:\n',path);
params
record.('Processing_Parameters') = params;

% Get list of images
[file_list, params] = update_list(target_dir, slash_str,params);

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

tag = num2str(round(datenum(datestr(clock))*100000));

cd(output_dir);
if isunix
    processed_str = ['/Foci_processed_', date,slash_str];
else
    processed_str = ['\Foci_processed_' , date, '_', tag,slash_str];
end


mkdir([output_dir processed_str]);
%tmp = [output_dir, '\',run_log_' date '_' name_out '.txt']
log_file = fopen([output_dir '\run_log_' date '_' name_out '.txt'], 'w');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%READ AND PROCESS:
num_file = length(file_list);
num_ch = size(well_struct.Label,2);

i = 1;
image_flag = 1;
while image_flag>0
    if mod(i,50)==0
        fprintf('Processed %d out of %d\n',i,num_file);
    end
    % Get updated list of images and modifiy num_file
    [file_list, params] = update_list(target_dir, slash_str,params);
    num_file = length(file_list);
    % THIS IF STATEMENT SHOULD BE CHANGED TO ADD A CHECKED FILE IN
    % DIRECTORY TELLING SOFTWARE IF IMAGING IS COMPLETE
    if i == num_file % if last image loaded, get out at end of loop
        image_flag=0;
    end
    try % This will write a log if image failed to process
        %% READING IMAGE
        nuc = readim([target_dir file_list{i}]);
        try
            for i_ch = 1:num_ch
                ch_img{i_ch} = readim([target_dir file_list{i}(1:end-8) well_struct.Label{i_ch} '.ics']);
            end
        catch
            % the file is not found - continue
            i = i+1;
            continue;
        end
        
        % This next line could be change by nuc_segmentor_local if we want
        % to optimize nuclear segmentation after acquisition. For now, it's
        % using the segmentation found during acquisition.
        mask = readim([target_dir file_list{i}(1:end-8) 'MASK.tif']);
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
            fprintf('segmenting nuclei in %s\n',[target_dir file_list{i}]);
        end
        fprintf(log_file,'segmenting nuclei in %s\n',[target_dir file_list{i}]);
        mask = erosion(dip_image(int16(mask)),1); % added on new pc on 7/7/2014
        if dimen >2
            mask = repmat(mask,1,1,size(ch_img{1},3));
            nuc = repmat(nuc,1,1,size(ch_img{1},3));
        end
        
        %%PROCESS SPOT DETECTION or Channel measurement
        
        %Initial nuc measurements
        try
            clear ms_nuc;
        end
        
        % measure general nuclear info (shape, size) and DAPI
        ms_nuc.image = repmat({file_list{i}(1:end-4)},1,num_nuc);
        if dimen>2
            ms_nuc_dapi = measure(squeeze(mask(:,:,0)),squeeze(nuc(:,:,0)),{'mean','size','p2a'}); % compute mean dapi for each nucleus
        else
            ms_nuc_dapi = measure(mask,nuc,{'mean','size','p2a'}); % compute mean dapi for each nucleus
        end
        ms_nuc.nuc_dapi = ms_nuc_dapi.mean;
        ms_nuc.nuc_area = ms_nuc_dapi.size*pix_size^2;
        ms_nuc.p2a = ms_nuc_dapi.p2a;
        
        % initialize record struct for db
        rec = [];
        rec.('nuc_dapi') = ms_nuc.nuc_dapi;
        rec.('nuc_area') = ms_nuc.nuc_area;
        rec.('p2a') = ms_nuc.p2a;
        rec.('nuc_dapi') = ms_nuc.nuc_dapi;
        image2analyze = [];
        image2analyze.('image_name') = file_list{i};
        image2analyze.('nuclei') = [];
        col_mask = nuc>-1; % initialize coloc mask to 1
        col_count = 0; % Will count the nubmer of time we overlap spot_mask
        % If col_count >1, then save coloc image
        
        %Start loop on different channels. Compute basic info (sum and
        %mean intensity) and if antibody is a spot antibody, compute spot
        %info
        for i_ch = 1:num_ch
            lab_str = well_struct.Antibody{i_ch};
            setfield(ms_nuc,sprintf('%s.background',lab_str),mean(ch_img{i_ch}(mask==0))); % indicator of background outside nucleus
            rec.(sprintf('%s.background',lab_str)) = mean(ch_img{i_ch}(mask==0));
            if dimen>2
                ms_nuc_ch = measure(squeeze(mask(:,:,0)),squeeze(mean(ch_img{i_ch},[],3)),{'mean','sum'});
            else
                ms_nuc_ch = measure(mask,ch_img{i_ch},{'mean','sum'});
            end
            setfield(ms_nuc,sprintf('%s.mean',lab_str),ms_nuc_ch.mean);
            setfield(ms_nuc,sprintf('%s.sum',lab_str),ms_nuc_ch.sum);
            rec.(sprintf('%s.mean',lab_str)) = ms_nuc_ch.mean;
            rec.(sprintf('%s.sum',lab_str)) = ms_nuc_ch.sum;
            %Check if antibody has spots. If yes, then enter spot
            %detection procedure
            if strcmp(lower(well_struct.Antibody{i_ch}),'53bp1') || strcmp(lower(well_struct.Antibody{i_ch}),'gh2ax')
                if PRINT_FLAG
                    fprintf('Spot detection\n');
                end
                [spot_struct,spot_mask,label_nuc,full_spot]= james_spot_detection7(ch_img{i_ch}, mask, params.max_foci_size, 100, params.min, params.k_val, 2, 1);
                col_mask = and(spot_mask,col_mask);
                col_count = col_count + 1;
                if sum(((spot_mask>0)*mask))>0 % make sure at least one focus was detected
                    % label all foci within nucleus by the nuc #
                    nuc_label_foci = dip_image(int16((spot_mask>0)*mask));
                    ms_nuc_foci = measure(nuc_label_foci,ch_img{i_ch},{'mean','size'}); % compute mean intensity of spot in its channel, and # of foci
                    setfield(ms_nuc,sprintf('%s.nfoci',lab_str),ms_nuc_foci.size);
                    setfield(ms_nuc,sprintf('%s.spot_mean',lab_str),ms_nuc_foci.mean);
                    rec.(sprintf('%s.nfoci',lab_str)) = ms_nuc_foci.size;
                    rec.(sprintf('%s.spot_mean',lab_str)) = ms_nuc_foci.mean;
                end
            end
            image2analyze.('nuclei') = [image2analyze.('nuclei'),rec];
            record.('images_analyzed') = [record.('images_analyzed'), image2analyze];
            
            % WRITE OUTPUT IMAGES
            %
            % Prepare output image
            % Project 3D images on 2D plane
            if size(nuc,3) > 1
                spot2D = squeeze(max(ch_img{i_ch},[],3));
                mask2D = squeeze(max(mask,[],3));
                nuc2D = squeeze(max(nuc,[],3));
                spot_mask2D = squeeze(max(spot_mask,[],3));
            else
                spot2D = ch_img{i_ch};
                mask2D = mask;
                nuc2D = nuc;
                spot_mask2D = spot_mask;
            end
            
            if strcmp(lower(well_struct.Antibody{i_ch}),'53bp1') || strcmp(lower(well_struct.Antibody{i_ch}),'gh2ax')
                if size(nuc,3) > 1
                    spot_mask2D = squeeze(max(spot_mask,[],3));
                else
                    spot_mask2D = spot_mask;
                end
            end
            %CREATE SPOT 2D DISPLAY
            spot2D = stretch(clip(spot2D,spot_max,spot_min)-spot_min); % Stretch all fitc images the same way
            overlay_img = overlay_phase(spot2D,stretch(nuc2D),'blue');
            overlay_img = overlay(overlay_img,xor(bdilation(spot_mask2D>0,2),bdilation(spot_mask2D>0)),[255 0 0]);
            edge = norm(gradient(mask2D))>0;
            edge = bskeleton(edge);
            overlay_img = overlay(overlay_img,edge,[255 255 0]);
            if print_label
                % Label nuclei with their numbers
                % Look up coordinates of nuc
                overlay_img
                msn = measure(mask2D,[],{'Minimum'});
                for i_label=1:size(msn,1)
                    try
                        text(max(msn(i_label).Minimum(1)-6,0),max(msn(i_label).Minimum(2)-6,0),sprintf('%d',i_label),'color','yellow','FontSize',10)
                    catch
                        fprintf('Nucleus %d overlap with another nucleus and it is lost in projection...\n',i_label);
                        fprintf(log_file,'Nucleus %d overlap with another nucleus and it is lost in projection...\n',i_label);
                    end
                end
                diptruesize(gcf,disp_size)
                saveas(gcf,[output_dir processed_str file_list{i}(1:end-4) '_' well_struct.Antibody{i_ch} '_spot_check.jpg'],'jpg');
                delete(gcf);
            else
                writeim(overlay_img,[output_dir processed_str file_list{i}(1:end-4) '_' well_struct.Antibody{i_ch} '_spot_check.jpg'],'jpeg');
            end
            
            %CREATE SPOT 2D DISPLAY WITHOUT FOCI CONTOUR
            overlay_img = overlay_phase(stretch(spot2D),stretch(nuc2D),'blue');
            overlay_img = overlay(overlay_img,edge,[255 255 0]);
            if print_label
                overlay_img
                % Label nuclei with their numbers
                for i_label=1:size(msn,1)
                    try
                        text(max(msn(i_label).Minimum(1)-6,0),max(msn(i_label).Minimum(2)-6,0),sprintf('%d',i_label),'color','yellow','FontSize',10)
                    catch
                        fprintf('Nucleus %d overlap with another nucleus and it is lost in projection...\n',i_label);
                        fprintf(log_file,'Nucleus %d overlap with another nucleus and it is lost in projection...\n',i_label);
                    end
                end
                diptruesize(gcf,disp_size)
                saveas(gcf,[output_dir processed_str file_list{i}(1:end-4) '_' well_struct.Antibody{i_ch} '_spot_check_no_foci_contour.jpg'],'jpg');
                delete(gcf);
            else
                writeim(overlay_img,[output_dir processed_str file_list{i}(1:end-4) '_' well_struct.Antibody{i_ch} '_spot_check_no_foci_contour.jpg'],'jpeg');
            end
        end
        if col_count > 1 % save coloc image and stats
            spot_mask2D = max(col_mask,[],3);
            overlay_img = overlay(overlay_img,xor(bdilation(spot_mask2D>0,2),bdilation(spot_mask2D>0)),[255 0 0]);
            edge = norm(gradient(mask2D))>0;
            edge = bskeleton(edge);
            overlay_img = overlay(overlay_img,edge,[255 255 0]);
            writeim(overlay_img,[output_dir processed_str file_list{i}(1:end-4) '_coloc_spot_check.jpg'],'jpeg');
            ms_nuc.col_nfoci = max(label(col_mask));
            rec.col_nfoci = ms_nuc.col_nfoci;
        end
        % WRITE OUTPUT FILE
        if (flag_first == 1) % print header with output if first time run
            nuc_sum_fh = fopen([output_dir processed_str 'nuc_summary_' name_out '.txt'], 'w');
            struct_print(nuc_sum_fh,ms_nuc,1);
            flag_first = 0;
        else
            nuc_sum_fh = fopen([output_dir processed_str 'nuc_summary_' name_out '.txt'], 'at');
            struct_print(nuc_sum_fh,ms_nuc,0);
        end
        fclose(nuc_sum_fh);
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

function [file_list, params] = update_list(target_dir, slash_str,params)

%get rid of directories in listing
file_list = dir([target_dir slash_str '*DAPI.tif']); %DAPI listed as tif on ASI
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
