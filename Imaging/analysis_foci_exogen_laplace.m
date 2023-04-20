function analysis_foci_exogen_laplace(target_dir, single_nuc_flag, params, PRINT_FLAG )
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

% for laplace, use 2D projection for detection. 3D is not so good (edge
% effects). Set dimen to 2.

if ~exist('PRINT_FLAG') || isempty(PRINT_FLAG); PRINT_FLAG = false; end

print_label = 0; % make it positive if you want to save images with label #
% Printing label requires to display each image on the
% screen (slow)
output_dir = target_dir;
if ~exist('name_out')
    name_out = ['foci_analysis_' date];
end
if ~exist('disp_size')
    disp_size = 100; % Set it to 100 if you want full image thumbnail. 50 if you want 50% (recommended)
end
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
dimen = 2;
flag_first = 1; % Flag used so that first image ready to have output print headers. Then set it to 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ANALYSIS:

tag = num2str(round(datenum(datestr(clock))*100000));

cd(output_dir);
if isunix
    processed_str = ['/Foci_processed_',date,'_', name_out,slash_str];
else
    processed_str = ['\Foci_processed_',date,'_', tag, '_', name_out,slash_str];
end


mkdir([output_dir processed_str]);
%tmp = [output_dir, '\',run_log_' date '_' name_out '.txt']
log_file = fopen([output_dir '\run_log_' date '_' name_out '.txt'], 'w');
param_file = fopen([output_dir '\param_used_' date '_' name_out '.txt'], 'w');
struct_print(param_file,params,1);
fclose(param_file);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%READ AND PROCESS:
num_file = length(file_list);
i = 1;
image_flag = 1;
while image_flag>0
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
            fitc_img = readim([target_dir file_list{i}(1:end-8) 'FITC.ics']);
        catch
            % the file is not found - continue
            i = i+1;
            continue;
        end
        
        if dimen==2
            fitc_img = squeeze(max(fitc_img,[],3));
            img = fitc_img;
        else
            img = squeeze(mean(fitc_img,[],3)); % Use image projection to detect nucleus
        end
        
        img_lap_nuc=threshold(img,'background');
        
        % Get rid of edges
        img_lap_nuc(:,0:1) = 0;
        img_lap_nuc(0:1,:) = 0;
        img_lap_nuc(:,end-1:end) = 0;
        img_lap_nuc(end-1:end,:) = 0;
        
        % take the envelop of the mask
        img_lap_nuc = convhull(img_lap_nuc,1);
        
        if dimen>2
            img_lap_nuc = repmat(img_lap_nuc,1,1,size(fitc_img,3));
            nuc = repmat(nuc,1,1,size(fitc_img,3));
        end
        mask = label(img_lap_nuc,2,100,10000000); % get rid of small isolated spots for mask
        num_nuc  = 1; %only one nuc
        
        %% Foci detection
        
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
        ms_nuc.nuc_area = ms_nuc_dapi.size;
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
        img_lap_spots=laplace(spot_img,1);               %tuning
        spot_mask=label(and(mask>0,img_lap_spots<-80),2,3,10000); 
        % keep spherical spots
        ms = measure(spot_mask,[],'p2a');
        if ~isempty(ms) % only remove non spherical foci if there is a focus
            spot_mask = label(and(spot_mask>0,msr2obj(spot_mask,ms,'p2a')<1.7));
        else
            spot_mask = newim(size(spot_img));
        end
        if sum(((spot_mask>0)*mask))>0 % make sure at least one focus was detected somewhere to compute any distance.
            % Add image name and nucID to structure output for
            % reference and create final structure output
            % To expedite programming, assumign only one nucleus per image
            % Should be true for Exogen
            ms_nuc.nfoci = max(spot_mask);
            ms_nuc.spot_mean = mean(spot_img(spot_mask>0));
            ms_nuc.perc_foci = sum(spot_mask>0)/sum(mask>0);
            ms_nuc.perc_foci_int = sum(spot_img(spot_mask>0))/sum(spot_img(mask>0));
        else % Else for if statement for testing at least one focus was detected.
            num_foci = 0;
            ms_foci = [];
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
        spot2D = stretch(clip(spot2D,spot_max,spot_min)-spot_min); % Stretch all fitc images the same way
        overlay_img = overlay_phase(spot2D,stretch(nuc2D),'blue');
        overlay_img = overlay(overlay_img,xor(spot_mask2D>0,berosion(spot_mask2D>0)),[255 0 0]);
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
            saveas(gcf,[output_dir processed_str file_list{i}(1:end-4) '_spot_check.jpg'],'jpg');
            delete(gcf);
        else
            writeim(overlay_img,[output_dir processed_str file_list{i}(1:end-4) '_spot_check.jpg'],'jpeg');
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
            saveas(gcf,[output_dir processed_str file_list{i}(1:end-4) '_spot_check_no_foci_contour.jpg'],'jpg');
            delete(gcf);
        else
            writeim(overlay_img,[output_dir processed_str file_list{i}(1:end-4) '_spot_check_no_foci_contour.jpg'],'jpeg');
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
