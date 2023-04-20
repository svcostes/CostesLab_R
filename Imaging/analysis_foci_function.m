function analysis_foci_function(target_dir,name_out,disp_size, PRINT_FLAG)
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

if ~exist('PRINT_FLAG') || isempty(PRINT_FLAG); PRINT_FLAG = false; end

param_dir = target_dir; %Directory where test_params.txt parameter file is located
output_dir = target_dir;
if ~exist('name_out')
    name_out = ['foci_analysis_' date];
end
if ~exist('disp_size')
    disp_size = 90; % Set it to 100 if you want full image thumbnail. 50 if you want 50% (recommended)
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
    target_dir = [target_dir,'\'];
end

% Error_list and batch_i are used in batch mode. If not set in case script
% is ran directly here, set them to 1 and load parameters
if ~exist('error_list','var')
    error_list = 1;
    batch_i = 1;
    %load parameters from the template filled out by user if not in batch
    %mode...
    params = load_params([param_dir 'test_params.txt']);
    %params = load_params('Z:\Hunter\decon\DNA_damage_staining\KO_gH2AX_Rad51\test_params.txt')
    %params.exp_str = exp_name;
else % If here, most likely it means script was ran already once or we are in batch mode
    if ~exist('keep_parameter') % if keep_parameter is not set, then we are not in batch mode, so load parameters
        params = load_params([param_dir 'test_params.txt']);
    else
        if keep_parameter == 0 % if keep_parameter exists, then we are in batch mode. Only load parameters if set that way
            params = load_params([param_dir 'test_params.txt']);
            %params.fitc_max_foci_size = 1;
        end
    end
end
% show parameters...
fprintf('Parameters settings:\n');
params
%buffer, from deconvolution:
buffer = 7;

% Set some parameters constant
if ~isfield(params,'fitc_max_foci_size')
    params.fitc_max_foci_size = 1;
end
if ~isfield(params,'txred_max_foci_size')
    params.txred_max_foci_size = 1;
end
params.target_dir = target_dir;
params.param_dir = param_dir;
params.output_dir = output_dir;
flag_first = 1; % Flag used so that first image ready to have output print headers. Then set it to 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ANALYSIS:
% Get file listing for directory. Depending on imaging platform
switch lower(params.media)
    case {'delta'}
        file_list = dir([target_dir '*.dv']);
    case {'meta'}
        file_list = dir([target_dir '*_D*.*']);
    case {'lsm'}
        file_list = dir([target_dir '*.lsm']);
    case {'zvi'}
        file_list = dir([target_dir '*.zvi']);
    case {'volocity'}
        file_list = dir([target_dir '*.ics']);
    case {'ics'}
        file_list = dir([target_dir '*.ics']);
    case {'tif'}
        file_list = dir([target_dir '*.tif']);
    otherwise
        file_list = dir(target_dir);
end


%get rid of directories in listing
ix = [file_list.isdir];
file_list = file_list(~ix);
file_list = {file_list.name};

%for output
dot_index = findstr(file_list{1}, '.');
params.suffix = file_list{1}(dot_index+1:end);
params.range_low = 1;
params.range_high = length(file_list);

cd(output_dir);
if isunix
    processed_str = ['/Foci_processed_',date,'_', name_out,slash_str];
else
    processed_str = ['\Foci_processed_',date,'_', name_out,slash_str];
end

mkdir([output_dir processed_str]);
log_file = fopen([output_dir 'run_log_' date '_' name_out '.txt'], 'w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%READ AND PROCESS:

for i = error_list(batch_i):length(file_list)
    %% READING IMAGE
    try
        switch lower(params.media)
            case {'delta'}
                [img,metadata] = bfopen3([target_dir file_list{i}]);
                % Make sure x,y scale not overkill. Anything below 0.15
                % um is set to 0.13
                if ~(isfield(params, 'dxy'))
                    params.dxy = metadata.get('Scale Factor for X');
                    params.dz = metadata.get('Scale Factor for Z');
                    if isempty(params.dxy) % For different deltavision version
                        params.dxy = metadata.get('Y element length (in um)');
                        params.dz = metadata.get('Z element length (in um)');
                    end
                    if ~(isfield(params, 'dxy'))
                        params.dxy = 1;
                    else
                        if params.dxy ==0
                            params.dxy = 1;
                        end
                    end
                    if ~(isfield(params, 'dz'))
                        params.dz = 1;
                    else
                        if params.dz ==0
                            params.dz = 1;
                        end
                    end
                end
                if params.dxy <0.11
                    xscale = 0.1/params.dxy;
                else
                    xscale = 1;
                end
                pix_size = max(0.1,params.dxy);
                zscale = pix_size/params.dz; % Always with respect to z since z scale always lower than xy
                % Make sure x,y scale not overkill. Anything below 0.15
                % um is set to 0.1
                nuc = img{params.dapi_index}(buffer:end-buffer,buffer:end-buffer,:);
                [w,h,d] = size(nuc);
                wsub = round(w/xscale); hsub = round(h/xscale);dsub = round(d/zscale);
                nuc = oversample(nuc,[wsub,hsub,dsub],'linear');
                nuc = squeeze(round(nuc(:,:,1:end))); % need this because oversample screws up first slice
                if (isfield(params, 'fitc_index'))
                    %fitc data is input. analyze this channel
                    fitc_img = img{params.fitc_index}(buffer:end-buffer,buffer:end-buffer,:);
                    fitc_img = oversample(fitc_img,[wsub,hsub,dsub],'linear');
                    fitc_img = squeeze(round(fitc_img(:,:,1:end))); % need this because oversample screws up first slice
                end
                if (isfield(params, 'txred_index'))
                    if ~isempty(params.txred_index)
                        %txred data is input. analyze this channel
                        txred_img = img{params.txred_index}(buffer:end-buffer,buffer:end-buffer,:);
                        txred_img = oversample(txred_img,[wsub,hsub,dsub],'linear');
                        txred_img = squeeze(round(txred_img(:,:,1:end)));% need this because oversample screws up first slice
                    end
                end
            case {'meta'}
                nuc = read_gray_stack([target_dir file_list{i}]);
                try
                    fitc_img = read_gray_stack([target_dir file_list{i}(1:end-8) 'F' file_list{i}(end-6:end)]);
                catch
                    %params=rmfield(params,'fitc_index');
                end
                try
                    txred_img = read_gray_stack([target_dir file_list{i}(1:end-8) 'T' file_list{i}(end-6:end)]);
                catch
                    %params=rmfield(params,'txred_index');
                end
            case {'volocity'}
                [img, metadata] = bfopen3([target_dir file_list{i}]);
                [bits, params.dxy, yscale, params.dz, channels, timeInt] = strread(metadata.get('parameter scale'));
                if ~(isfield(params, 'dxy'))
                    params.dxy = 1;
                else
                    if params.dxy ==0
                        params.dxy = 1;
                    end
                end
                if ~(isfield(params, 'dz'))
                    params.dz = 1;
                else
                    if params.dz ==0
                        params.dz = 1;
                    end
                end
                % Make sure x,y scale not overkill. Anything below 0.15
                % um is set to 0.13
                if params.dxy <0.11
                    xscale = 0.1/params.dxy;
                else
                    xscale = 1;
                end
                pix_size = max(0.1,params.dxy);
                zscale = pix_size/params.dz; % Always with respect to z since z scale always lower than xy
                nuc = gaussf(img{params.dapi_index}(buffer:end-buffer,buffer:end-buffer,:)); %this will smooth these images. Too pixelated...
                [w,h,d] = size(nuc);
                wsub = round(w/xscale); hsub = round(h/xscale);dsub = round(d/zscale);
                nuc = oversample(nuc,[wsub,hsub,dsub],'linear');
                nuc = round(nuc(:,:,1:end)); % need this because oversample screws up first slice
                if (isfield(params, 'fitc_index'))
                    %fitc data is input. analyze this channel
                    fitc_img = gaussf(img{params.fitc_index}(buffer:end-buffer,buffer:end-buffer,:));
                    fitc_img = oversample(fitc_img,[wsub,hsub,dsub],'linear');
                    fitc_img = round(fitc_img(:,:,1:end)); % need this because oversample screws up first slice
                end
                if (isfield(params, 'txred_index'))
                    if ~isempty(params.txred_index)
                        %txred data is input. analyze this channel
                        txred_img = gaussf(img{params.txred_index}(buffer:end-buffer,buffer:end-buffer,:));
                        txred_img = oversample(txred_img,[wsub,hsub,dsub],'linear');
                        txred_img = round(txred_img(:,:,1:end));% need this because oversample screws up first slice
                    end
                end
            case {'tif'}
                img = readim([target_dir file_list{i}]);
                nuc = img{params.dapi_index};
                if (isfield(params, 'fitc_index'))
                    %fitc data is input. analyze this channel
                    fitc_img = img{params.fitc_index};
                end
                if (isfield(params, 'txred_index'))
                    %txred data is input. analyze this channel
                    txred_img = img{params.txred_index};
                end
                pix_size = 0.1;
            case {'ics'}
                [img, metadata] = bfopen3([target_dir file_list{i}]);
                if length(size(img))==2 % Images is a tensor. For ics, typically 2D color
                    nuc = img{params.dapi_index};
                    if (isfield(params, 'fitc_index'))
                        %fitc data is input. analyze this channel
                        fitc_img = img{params.fitc_index};
                    end
                    if (isfield(params, 'txred_index'))
                        %txred data is input. analyze this channel
                        txred_img = img{params.txred_index};
                    end
                else
                    if length(size(img))==3 % probably 2D image
                        img = reshape(img,[size(img,1),size(img,2),1,size(img,3)]);
                    end
                    nuc = squeeze(img(:,:,:,params.dapi_index-1));
                    if (isfield(params, 'fitc_index'))
                        %fitc data is input. analyze this channel
                        fitc_img = squeeze(img(:,:,:,params.fitc_index-1));
                    end
                    if (isfield(params, 'txred_index'))
                        %txred data is input. analyze this channel
                        txred_img = squeeze(img(:,:,:,params.txred_index-1));
                    end
                end
                if ~(isfield(params, 'dxy'))
                    params.dxy = 1;
                else
                    if params.dxy ==0
                        params.dxy = 1;
                    end
                end
                if ~(isfield(params, 'dz'))
                    params.dz = 1;
                else
                    if params.dz ==0
                        params.dz = 1;
                    end
                end
                % Make sure x,y scale not overkill. Anything below 0.15
                % um is set to 0.13
                if params.dxy <0.11
                    xscale = 0.1/params.dxy;
                else
                    xscale = 1;
                end
                pix_size = max(0.1,params.dxy);
                zscale = pix_size/params.dz; % Always with respect to z since z scale always lower than xy
                
            otherwise
                %%% assume otherwise case = zvi file!!!
                [img, metadata] = bfopen3([target_dir file_list{i}]);
                nuc = img{params.dapi_index};
                if (isfield(params, 'fitc_index'))
                    %fitc data is input. analyze this channel
                    fitc_img = img{params.fitc_index};
                end
                if (isfield(params, 'txred_index'))
                    %txred data is input. analyze this channel
                    txred_img = img{params.txred_index};
                end
                if ~(isfield(params, 'dxy'))
                    params.dxy = 1;
                else
                    if params.dxy ==0
                        params.dxy = 1;
                    end
                end
                if ~(isfield(params, 'dz'))
                    params.dz = 1;
                else
                    if params.dz ==0
                        params.dz = 1;
                    end
                end
                % Make sure x,y scale not overkill. Anything below 0.15
                % um is set to 0.13
                if params.dxy <0.11
                    xscale = 0.1/params.dxy;
                else
                    xscale = 1;
                end
                pix_size = max(0.1,params.dxy);
                zscale = pix_size/params.dz; % Always with respect to z since z scale always lower than xy
        end
        if (length(size(nuc)) == 3)
            dimen = 3;
        else
            dimen = 2;
        end
        
        %% FIX ILLUMINATION PROBLEM
        %         if (dimen == 3)
        %             % Fix illumination problem
        %             tmp = gaussf(max(nuc,[],3),200);
        %             nuc = nuc - repmat(squeeze(tmp),1,1,size(nuc,3));
        %         else
        %             tmp = gaussf(nuc,200);
        %             nuc = nuc - tmp;
        %         end
        %%SEGMENT NUCLEI
        if PRINT_FLAG
            fprintf('dxy=%f,dz=%f,xscale=%f,zscale=%f\n',params.dxy,params.dz,xscale,zscale);
            fprintf('segmenting nuclei in %s\n',[target_dir file_list{i}]);
        end
        fprintf(log_file,'segmenting nuclei in %s\n',[target_dir file_list{i}]);
        if isempty(params.nuc_th)
            th_level = 1.5;
        else
            th_level = params.nuc_th;
        end
        if dimen==3
            mask = nuc_segmentor_local(max(nuc,[],3),params.nuc_rad,1,th_level,1,0); %threshold level is 1.7
        else
            mask = nuc_segmentor_local(nuc,params.nuc_rad,1,th_level,1,0); %threshold level is 1.7
        end
        % get rid of small nuclei or deformed nuclei (values
        % predetermined for ImageXpress images)
        ms = measure(mask,[],{'size','p2a'});
        mask_p2a = msr2obj(mask,ms,'p2a');
        mask_size = msr2obj(mask,ms,'size');
        norm_size = pi()*(params.nuc_rad*1.5)^2; %needs 1.5 because segmentation bias towards larger size
        keep_mask = and(and(mask_p2a<1.35,mask>0),and(mask_size>norm_size*0.5,mask_size<norm_size*1.75));
        mask = mask*keep_mask;
        %         mask = dilation(erosion(mask,params.nuc_rad*1.3),params.nuc_rad*1.3);
        mask = relabel_v2(mask);
        if dimen == 3
            mask = repmat(mask,1,1,size(nuc,3));
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
        
        num_nuc = max(mask);
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
        if sum(spot_mask)>0 % make sure at least one focus was detected somewhere to compute any distance.
            % label all foci within nucleus by the nuc #
            nuc_label_foci = dip_image(uint8((spot_mask>0)*mask));
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
            nuc_sum_fh = fopen([output_dir processed_str 'nuc_summary_' name_out '.txt'], 'w');
            foci_sum_fh = fopen([output_dir processed_str 'foci_summary_' name_out '.txt'], 'w');
            struct_print(nuc_sum_fh,ms_nuc,1);
            struct_print(foci_sum_fh,ms_foci,1);
            flag_first = 0;
        else
            nuc_sum_fh = fopen([output_dir processed_str 'nuc_summary_' name_out '.txt'], 'at');
            foci_sum_fh = fopen([output_dir processed_str 'foci_summary_' name_out '.txt'], 'at');
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
            %             fitc2D = squeeze(max(fitc_img,[],3));
            %             txred2D = squeeze(max(txred_img,[],3));
            spot2D = squeeze(max(spot_img,[],3));
            mask2D = squeeze(max(mask,[],3));
            nuc2D = squeeze(max(nuc,[],3));
            spot_mask2D = squeeze(max(spot_mask,[],3));
        else
            %             fitc2D = fitc_img;
            %             txred2D = txred_img;
            spot2D = spot_img;
            mask2D = mask;
            nuc2D = nuc;
            spot_mask2D = spot_mask;
        end
        % Look up coordinates of nuc
        msn = measure(mask2D,[],{'Minimum'});
        
        %CREATE SPOT 2D DISPLAY
        spot2D = stretch(clip(spot2D,spot_max,spot_min)-spot_min); % Stretch all fitc images the same way
        overlay_img = overlay_phase(stretch(spot2D),stretch(nuc2D),'blue');
        overlay_img = overlay(overlay_img,xor(bdilation(spot_mask2D>0,2),bdilation(spot_mask2D>0)),[255 0 0]);
        edge = norm(gradient(mask2D))>0;
        edge = bskeleton(edge);
        overlay_img = overlay(overlay_img,edge,[255 255 0])
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
        saveas(gcf,[output_dir processed_str file_list{i}(1:end-4) '_spot_check.jpg'],'jpg');
        delete(gcf);
        
        %CREATE SPOT 2D DISPLAY WITHOUT FOCI CONTOUR
        overlay_img = overlay_phase(stretch(spot2D),stretch(nuc2D),'blue');
        overlay_img = overlay(overlay_img,edge,[255 255 0])
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
        
        %CREATE SPOT 2D DISPLAY WITHOUT FOCI CONTOUR and getting rid of
        %background (enhanced spot image)
        if size(nuc,3) > 1
            spot2D = squeeze((mask(:,:,0)>0)*(max(spot_img,[],3)-gaussf(mean(spot_img,[],3),7)));
        else
            spot2D = clip(spot_img - gaussf(spot_img,7),1e5,0);
        end
        
        overlay_img = overlay_phase(stretch(spot2D),stretch(nuc2D),'blue')
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
        saveas(gcf,[output_dir processed_str file_list{i}(1:end-4) '_enhanced_spot_nobackground.jpg'],'jpg');
        delete(gcf);
        
        
    catch % Try to process image. If failed, write log.
        error_struct = lasterror;
        disp(['ERROR! ' error_struct.identifier '\n']);
        disp([error_struct.message '\n']);
        disp('Traceback:')
        fprintf(log_file,'ERROR! %s\n%s\nTraceback:\n', error_struct.identifier, error_struct.message);
        for err_i=1:length(error_struct.stack)
            disp(['name: ' error_struct.stack(err_i).name ', line: ' num2str(error_struct.stack(err_i).line)]);
            fprintf(log_file, 'name: %s, line: %d\n', error_struct.stack(err_i).name, error_struct.stack(err_i).line);
        end
    end
end

fclose(log_file);