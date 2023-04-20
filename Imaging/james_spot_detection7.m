function [spot_struct,spot_mask,label_nuc,full_spot,trim_full_spot,sphase_indices, return_wps_arr] = james_spot_detection7(spot, label_nuc, max_size, max_int, min_int, k_val, res, gen_wps_flag, wps_arr, nuc_coords, PRINT_FLAG)
%
% [spot_struct,spot_mask,label_nuc,full_spot,trim_full_spot,sphase_indices,
% return_wps_arr] =
% james_spot_detection6_sc(spot, label_nuc, max_size, max_int, min_int,
% k_val, res, gen_wps_flag, wps_arr, nuc_coords)
%
% This function will take a nuclear segmentation and it's corresponding
% spot signal and detect and quantify each spot information
%
% Input:
%       spot: gray image of spot
%       label_nuc: labeled image (segmentation)
%       max_size: Size used for dip_localminima algorithm
%       max_int: Intensity used for dip_localminima. Not very sensitive. Leave blank (i.e. [] to use default value: 100).
%       min_int: Absolute intensity cutoff for foci on maximum intensity in full focus. Anythig below is
%       rejected at the end of the detection process.
%       k_val, res: wavelets parameters for spot detection
%       gen_wps_flag: flag for wavelet. Set to 1 by default. If set to 0,
%       then need to enter previous computed wavelet (to save time). This
%       is never used. 
%       wps_arr: only if gen_wps_flag set to 0.
%       nuc_coords: (optional) Coordinates delimiting the nucleus (used
%       when image was previously cropped, to save time)

%
%
% Output:
%       spot_struct:
%           'nucID',1:max(label_nuc),...
%           'nuc_coords' % coords of center of nucleus
%           'count' % foci/nucleus
%           'foci' % [nucID, fociID, foci_size,foci_meanI,foci_maxI];
%       spot_label:labeled center foci
%       label_nuc: labeled nuclei
%       full_spot: labeled full foci
%       trim_full_spot: labeled full foci, only showing in their center
%       slice (smaller)  *** MAKE IT OPTIONAL
%
%
% Copyright: Lawrence Berkeley National Laboratory
% Author: Sylvain Costes & James Chen
% November 2005
% Modified November 2009
% Modified Aug 2010
% Modified March 2014, by Sylvain. Changed min_int cutoff on foci to use
% the maximum intensity found in the spot as needing to be above the
% cutoff.

if ~exist('PRINT_FLAG') || isempty(PRINT_FLAG); PRINT_FLAG = false; end

dimen = length(size(spot));     % Determine if dimension of image (i.e. 3D or 2D)

if ~exist('max_size')           %avg size of foci set to 40 or 60 for 2D or 3D pictures,
    max_size = 5;           %default value is 5
end

if ~exist('max_int')            %avg intensity set to 1400 or 2100 for 2D or 3D pictures,
    max_int = 100;        %default value is 100
end

if dimen == 3                   %if stack, get number of slices
    nslice = size(spot,3);
    nslice2 = size(label_nuc,3);
    if nslice~=nslice2
        error('spot and label_nuc images do not have the same dimensions');
    end
end

min_size = 10; % empirical value to delete final foci that are too small
%--------------------------------------------------------------------------

%ID gaps were closed in nuc_segmentor
total_nuc = max(label_nuc); %number of nuclear areas isolated

% Get coordinates of all nuclei.
keep_nuc = [];
info_spot = [];
wave_foci_mask = newim(size(spot));


% try
%     if size(nuc_coords,1) ~= total_nuc
%         error('The labelled nuc image must have the same number of nuclei as the coordinates entered');
%     end
%     keep_nuc = 1:total_nuc;
% catch %If gets here, means that nuc_coords was not entered. Then compute it on a subs
%     fprintf('Initial # of nuclei:%d\n',total_nuc);
%
%     if dimen == 3
%         %find which z slice has the greatest area for this nuc id
%         %first find out the range of z slices that this particular nuc id spans
%         temp_sum_label_nuc = sum(label_nuc, [], 2);
%         temp_sum_label_nuc = sum(temp_sum_label_nuc, [], 1);
%         [max_area, max_z] = max(temp_sum_label_nuc);
%         max_z = max_z(3);
%
%         %old one. keep in case.
%         %small_nuc = squeeze(label_nuc(0:2:end,0:2:end,round(end/2)));
%         small_nuc = squeeze(label_nuc(0:2:end,0:2:end,max_z));
%     else
%         small_nuc = label_nuc(0:2:end,0:2:end);
%     end
%     [small_w, small_h] = size(small_nuc);
%     %need to remove gaps in ids, lost from the subsampling the smaller nuc picture and re adjust total_nuc
%     h = waitbar(0, 'Small Nuc Adjustment...', 'position', [100 100 270 56.25]);
%     temp_counter = 1;
%     for p=1:max(small_nuc)
%         waitbar(p/max(small_nuc));
%         if ~isempty(find(small_nuc==p))
%             small_nuc(find(small_nuc==p)) = temp_counter;
%             temp_counter = temp_counter + 1;
%         end
%     end
%     close(h);
%
%     total_nuc = max(small_nuc);
%
%     h = waitbar(0,'Computing nuclear limits...', 'position', [100 100 270 56.25]);
%     for nuc_num = 1:total_nuc
%         waitbar(nuc_num/total_nuc);
%
%         [temp_img,coords] = crop_from_mask(small_nuc,small_nuc==nuc_num);
%         if isempty(coords)
%             fprintf('haha, empty!\n');
%         end
%         if coords(1) < 0
%             coords(1) = 0;
%         end
%         if coords(2) < 0
%             coords(2) = 0;
%         end
%         if coords(3) < 0
%             coords(3) = 0;
%         end
%         if coords(4) < 0
%             coords(4) = 0;
%         end
%
%         if dimen == 3
%             nuc_coords(nuc_num,:,:) = [2.*coords;[0 nslice-1]];
%         else
%             nuc_coords(nuc_num,:,:) = 2.*coords;
%         end
%         diff_coords = squeeze(abs(nuc_coords(nuc_num,:,1)-nuc_coords(nuc_num,:,2)));
%         cropped_area = diff_coords(1)*diff_coords(2);
%         keep_nuc = [keep_nuc,nuc_num];
%     end
%     close(h);
%     %need to remove gaps in the smaller nuc picture and re adjust total_nuc
%
%     nuc_coords = nuc_coords(keep_nuc,:,:);
%     total_nuc = length(keep_nuc);
%
%     if (dimen == 2)
%         label_nuc = squeeze(label_nuc);
%     end
%
%     fprintf('\nKept %d nuclei.\n\n',total_nuc);
% end

foci_mask = newim(size(spot));
trim_foci_mask = newim(size(spot));
spot_single = newim(size(spot));
sphase_indices = [];
return_wps_arr = [];
if PRINT_FLAG
    fprintf('NUC NUM: ');
end
for nuc_num = 1:total_nuc           %parse through all nuclei
    if PRINT_FLAG
        fprintf('%d ... ', nuc_num);
    end
    %    coords = squeeze(nuc_coords(nuc_num,:,:));          %grab coordinates of local nucleus we are analyzing
    if dimen == 3
        %% 3D
        [local_nuc,coords] = crop_from_mask(label_nuc,label_nuc==nuc_num);     %square crop around nuc in question
        
        % with new cropped implementation this can be empty
        if isempty(coords)
            continue;
        end
        
        nuc_coords(nuc_num,:,:) = coords;
        local_spot = spot(coords(1,1):coords(1,2),coords(2,1):coords(2,2),coords(3,1):coords(3,2));     %square crop around original spot image of nuc in question
        nuc_labels = unique(uint8(local_nuc));                          %get unique nuc ids within local_nuc (on edges of square crop
        nuc_labels_mask = newim(size(local_spot));                      %create new blank image with same size as local_nuc
        nuc_labels_mask(:,:,:) = 1;                                     %set all values to nuc_labels_mask to 1
        for y = 1:size(nuc_labels,1)                                    %for each unique nuc ids within local_nuc
            if (nuc_labels(y) ~= 0 && nuc_labels(y) ~= nuc_num)         %if the id is not 0 (background) or the nuc_num in question, set it to 0 in nuc_labels_mask
                nuc_labels_mask(find(local_nuc==nuc_labels(y))) = 0;
            end
        end

        % Clean up image and make it less intensity dependent by
        % subtracting blurred local_spot
        % Commented out on Aug 26, 2010
%         cor_spot = newim(size(local_spot));
%         for i=1:size(local_spot,3)
%             cor_spot(:,:,i-1) = gaussf(local_spot(:,:,i-1),1)-gaussf(local_spot(:,:,i-1),6);
%         end
        cor_spot = local_spot;
       
        local_nuc = local_nuc*(local_nuc==nuc_num);                     %isolate only the nuc_num id in local_nuc

                if gen_wps_flag
                    try
                        [local_foci_mask_r2, wps_r2] = wavelet_spot_detection(cor_spot, local_nuc>0, 2, (1.5*k_val));
                        [local_foci_mask_r3, wps_r3] = wavelet_spot_detection(cor_spot, local_nuc>0, 3, k_val);
                        local_foci_mask = or(local_foci_mask_r2, local_foci_mask_r3);
                        %[local_foci_mask, wps] = wavelet_spot_detection(cor_spot, local_nuc>0, res, k_val);
                    catch
                        fprintf('Error with wavelet\n');
                    end
                else
                    [local_foci_mask, wps] = wavelet_spot_detection(cor_spot, local_nuc>0, res, k_val, wps_arr{nuc_num});     %use wavelet algorithm to generate local mask of foci
                end
                
                % If image is saturated, saturated spots have a hole inside
                % the flat intensity part. Fill it back as it is an
                % artifact of wavelet
                local_foci_mask = (map('berosion',map('bdilation',local_foci_mask,'4')>0,'3'))>0;
                if (max(local_spot)<256)
                    local_foci_mask = or(local_foci_mask,local_spot==255);
                else
                    local_foci_mask = or(local_foci_mask,local_spot==4095);
                end
                local_foci_mask = and(local_foci_mask>0,local_nuc>0);  % make sure local_foci_mask is in local_nuc
                local_foci_mask = label(local_foci_mask>0,3,min_size,100000000)>0;
%         ms = measure(threshold(cor_spot),cor_spot,{'maxval','size'}); % measure properties of obvious foci (isodata)
%         ind = find(ms.size>5); % isolate foci that are real
%         local_foci_mask = cor_spot> 0.5*min(ms.maxval(ind)); % Set 1/2 of max of dimmest foci as cut off. Then resegment local spot.
 %        local_foci_mask = object_separation(local_foci_mask,1)>0; % Get rid of small object and separate some clusters

    %changed by james on 12/10/09. object_separation seems to be the reason
    %why bright foci are not being picked up. it is identified after the
    %wavelet detection, but object separation is clearing it out. gonna try
    %it with a gray scaled image passed in to get a more accurate
    %grow_region. UPDATE:
    
    %[local_foci_mask_os, local_spot_single] = object_separation_jc(local_foci_mask,0.1, cor_spot); % Get rid of small object and separate some clusters
    %local_foci_mask = local_foci_mask>0;
    
    %temp_holder = newim(size(wave_foci_mask));
        %temp_holder(coords(1,1):coords(1,2),coords(2,1):coords(2,2),coords(3,1):coords(3,2)) = local_foci_mask;
        %wave_foci = wave_foci + temp_holder;

        %james 12/10/09: this section is essentially done by
        %object_separation. we just need to save the spots, so i modified
        %object_separation to return the seed image too.
        %%%%%%%%%%%%%%%%%%%%%%commenting out, 12/10/09%%%%%%%%%%%%%%%%%%%%
        %for debugging:
        cor_spot = clip(cor_spot, 50000,0);
        temp_cor_spot = newim(size(cor_spot,1), size(cor_spot,2), size(cor_spot,3)+2);
        temp_cor_spot(:,:,1:end-1) = cor_spot;
        temp_local_foci_mask = bdilation(local_foci_mask);
        temp_local_foci_mask = newim(size(local_foci_mask,1), size(local_foci_mask,2), size(local_foci_mask,3)+2);
        temp_local_foci_mask(:,:,1:end-1) = local_foci_mask;
        
        
        if max(local_foci_mask > 0)        %if the local foci mask is not empty
            %apparently, dip_localminima only looks at z indices 1:end-1
            %local_spot_single = dip_localminima(-(add_z_ends(cor_spot)), add_z_ends(bdilation(local_foci_mask))>0, 3, max_int, max_size, 1);   %initial foci seed detection
            %local_spot_single = dip_localminima(-temp_cor_spot, temp_local_foci_mask>0, 3, max_int, max_size, 1);   %initial foci seed detection
            try
                local_spot_single = dip_localminima(-(add_z_ends(local_spot)), add_z_ends(local_foci_mask)>0, 3, max_int, max_size, 1);
                %local_spot_single = dip_localminima(-(add_z_ends(local_spot)), add_z_ends(local_foci_mask)>0, 3, max_int, max_size, 1);   %initial foci seed detection
                local_spot_single = local_spot_single(:,:,1:end-1);
                local_foci_mask = dip_growregions(label(local_spot_single),[],local_foci_mask, dimen, 0, 'high_first');
                % Weed out dim foci (mean Int < min_int)
                % S. Costes, September 2009
                % trying with mean Int < min_int + 0.1(max_int-min_int)
                temp_ms = measure(local_foci_mask, local_spot, {'MaxVal','size'}, [], 2);
                temp_mean = msr2obj(local_foci_mask,temp_ms,'MaxVal');
                temp_size = msr2obj(local_foci_mask,temp_ms,'size');
                local_foci_mask = and(temp_mean>min_int,temp_size>min_size);
            catch
                local_spot_single = newim(size(local_foci_mask));
            end
        else
            local_spot_single = newim(size(local_foci_mask));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %make it a binary        
        local_foci_mask = local_foci_mask>0;

        
        local_spot_single = (local_spot_single*local_nuc*local_foci_mask)>0;            %isolate seeds only within the nuclear segmented area/volume, also in the local_foci_mask
        local_label_spot_single = label(local_spot_single,3);
        %need to deal with grouped seeds
        local_label_spot_single_ms = measure(local_label_spot_single, [], {'size', 'center'}, [], 3);  %take measurement of size and center of labeled local_spot_single
        if isempty(local_label_spot_single_ms)
            large_seeds = [];
        else
            large_seeds = find(local_label_spot_single_ms.size>1);      %find spots > 1 in size
        end
        if ~isempty(large_seeds)        %if the large seed list isnt' empty
            for m = large_seeds         %for each id associated with a seed > 1 in size
                new_large_seed_coords = round(local_label_spot_single_ms(m).center);    %grab the center coordinates and round it
                local_label_spot_single(find(local_label_spot_single==m)) = 0;          %set the original pixels to 0
                local_label_spot_single(new_large_seed_coords(1), new_large_seed_coords(2), new_large_seed_coords(3)) = m;      %set the new rounded coordinate to the id
            end
        end


        local_spot_single = local_label_spot_single>0;      %the labeled version has been filtered for seeds > 1. take the binary and set it back to local_spot_single


        %grow regions back
        %local_label_foci_mask = dip_growregions(local_label_spot_single, [], local_foci_mask, 3, 0, 'high_first', 'default');
        local_label_foci_mask = dip_growregions(local_label_spot_single, [], local_foci_mask, 3, 0, 'high_first');

        %OK TO HERE%
        
        %%%added 12/20/07 - test gravity method of obtaining spot_mask.
        %%for each of the foci (in 3D), find the center slice (based off of
        %%ms.Center) and get rid of all the other slices for this foci.
        %%this "trimmed" down version of the local labeled foci mask is
        %%what is used.
        local_spot_single_ms = measure(dip_image(uint8(label(local_spot_single,1))), [], 'center', [], 2);     %get ids of all foci, and centers
        trim_local_label_foci_mask = newim(size(local_label_foci_mask));
        %[h w] = size(trim_local_label_spot_mask);
        if ~isempty(local_spot_single_ms)
            for k = 0:size(local_label_foci_mask,3)-1      %for each slice
                %fprintf('k: %d\n', k);

                indices = find(round(local_spot_single_ms.Center(3,:))==k);    %find the foci that have the z centers of the foci in this slice.
                keep_ids = local_spot_single_ms.ID(indices);
                if ~isempty(keep_ids)

                    %HERE NEED TO FIND A WAY TO WEED OUT BAD
                    %SLICES THAT DON'T WANT.
                    temp_llfm = double(local_label_foci_mask(:,:,k));       %isolate the slice, translate to double array
                    keep_coords = compare2Arrays(temp_llfm(:), keep_ids);   %get coordinates of the ids we want to keep.
                    temp_trim_local_label_foci_mask = newim(size(squeeze(local_label_foci_mask(:,:,k))));
                    temp_trim_local_label_foci_mask(keep_coords-1) = 1;     %create binary mask of the foci mask for this slice
                    temp_trim_local_label_foci_mask = temp_trim_local_label_foci_mask*squeeze(local_label_foci_mask(:,:,k));        %slice of k, that only has the masks of the foci with centers in this slice

                    %spot_center = measure(dip_image(uint8(temp_trim_local_label_foci_mask)), (squeeze(local_spot(:,:,k))), 'gravity', []);
                    %spot_center_coords = round(spot_center.Gravity)';
                    %trim_crop_fitc_spot(

                    %trim_local_label_spot_mask(spot_center_coords(:,1)*w+spot_center_coords(:,2)+k*w*h) = 1;

                    trim_local_label_foci_mask(:,:,k) = temp_trim_local_label_foci_mask;
                end

                %trim_crop_fitc_foci(:,:,k) = crop_fitc_foci(:,:,k) == keep_ids;
            end
        end
        
        %HERE WE HAVE THE LOCAL LABEL FOCI MASK, THE TRIMMED LOCAL LABEL
        %FOCI MASK AND THE SPOT SINGLE MASK. CHECK WHICH FOCI POSITIONS ARE
        %BETTER BY OVERLAYING WITH THE LOCAL SPOT IMAGE
        %%%



        %grab foci measurements
        if (max(local_label_foci_mask) == 0)
            count_spot(nuc_num) = 0;
            info_spot = [info_spot; nuc_num, NaN, NaN, NaN, NaN];
        else
            %take measurements for foci information
            foci_ms = measure(local_label_foci_mask,local_spot,{'Mean','MaxVal','Size'}, [], 3);
            spot_mean = foci_ms.mean;
            spot_max = foci_ms.maxval;
            spot_size = foci_ms.size;
            spot_ID = foci_ms.ID;
            num_spot = length(spot_mean);
            info_spot = [info_spot; ones(num_spot,1)*nuc_num, spot_ID', spot_size', spot_mean', spot_max'];
            count_spot(nuc_num) = num_spot;
        end

        %patch local_label_foci_mask into foci_mask
        foci_mask(coords(1,1):coords(1,2),coords(2,1):coords(2,2),coords(3,1):coords(3,2)) = foci_mask(coords(1,1):coords(1,2),coords(2,1):coords(2,2),coords(3,1):coords(3,2)) + local_label_foci_mask;

        %patch trim_local_label_foci_mask into trim_foci_mask
        trim_foci_mask(coords(1,1):coords(1,2),coords(2,1):coords(2,2),coords(3,1):coords(3,2)) = trim_foci_mask(coords(1,1):coords(1,2),coords(2,1):coords(2,2),coords(3,1):coords(3,2)) + trim_local_label_foci_mask;

        %patch local_spot_single into spot_single
        spot_single(coords(1,1):coords(1,2),coords(2,1):coords(2,2),coords(3,1):coords(3,2)) = spot_single(coords(1,1):coords(1,2),coords(2,1):coords(2,2),coords(3,1):coords(3,2)) + local_spot_single;


    else
        %% 2D
        [local_nuc,coords] = crop_from_mask(label_nuc,label_nuc==nuc_num);     %square crop around nuc in question
        nuc_coords(nuc_num,:,:) = coords;
        local_spot = spot(coords(1,1):coords(1,2),coords(2,1):coords(2,2));     %square crop around original spot image of nuc in question
        nuc_labels = unique(uint8(local_nuc));                          %get unique nuc ids within local_nuc (on edges of square crop
        nuc_labels_mask = newim(size(local_spot));                      %create new blank image with same size as local_nuc
        nuc_labels_mask(:,:) = 1;                                     %set all values to nuc_labels_mask to 1
        for y = 1:size(nuc_labels,1)                                    %for each unique nuc ids within local_nuc
            if (nuc_labels(y) ~= 0 && nuc_labels(y) ~= nuc_num)         %if the id is not 0 (background) or the nuc_num in question, set it to 0 in nuc_labels_mask
                nuc_labels_mask(find(local_nuc==nuc_labels(y))) = 0;
            end
        end

        % Clean up image and make it less intensity dependent by
        % subtracting blurred local_spot
        cor_spot = newim(size(local_spot));
        cor_spot = gaussf(local_spot,1)-gaussf(local_spot,6);

        local_nuc = local_nuc*(local_nuc==nuc_num);                     %isolate only the nuc_num id in local_nuc

                if gen_wps_flag
                    try
                        [local_foci_mask_r2, wps_r2] = wavelet_spot_detection(cor_spot, local_nuc>0, 2, (1.5*k_val));
                        [local_foci_mask_r3, wps_r3] = wavelet_spot_detection(cor_spot, local_nuc>0, 3, k_val);
                        %[local_foci_mask, wps] = wavelet_spot_detection(cor_spot, local_nuc>0, res, k_val);
                    catch
                        fprintf('Error with wavelet\n');
                    end
                else
                    [local_foci_mask, wps] = wavelet_spot_detection(cor_spot, local_nuc>0, res, k_val, wps_arr{nuc_num});     %use wavelet algorithm to generate local mask of foci
                end
                
                local_foci_mask = or(local_foci_mask_r2, local_foci_mask_r3);
%         ms = measure(threshold(cor_spot),cor_spot,{'maxval','size'}); % measure properties of obvious foci (isodata)
%         ind = find(ms.size>5); % isolate foci that are real
%         local_foci_mask = cor_spot> 0.5*min(ms.maxval(ind)); % Set 1/2 of max of dimmest foci as cut off. Then resegment local spot.
 %        local_foci_mask = object_separation(local_foci_mask,1)>0; % Get rid of small object and separate some clusters

    %changed by james on 12/10/09. object_separation seems to be the reason
    %why bright foci are not being picked up. it is identified after the
    %wavelet detection, but object separation is clearing it out. gonna try
    %it with a gray scaled image passed in to get a more accurate
    %grow_region. UPDATE:
    
    %[local_foci_mask_os, local_spot_single] = object_separation_jc(local_foci_mask,0.1, cor_spot); % Get rid of small object and separate some clusters
    %local_foci_mask = local_foci_mask>0;
    
    %temp_holder = newim(size(wave_foci_mask));
        %temp_holder(coords(1,1):coords(1,2),coords(2,1):coords(2,2),coords(3,1):coords(3,2)) = local_foci_mask;
        %wave_foci = wave_foci + temp_holder;

        %james 12/10/09: this section is essentially done by
        %object_separation. we just need to save the spots, so i modified
        %object_separation to return the seed image too.
        %%%%%%%%%%%%%%%%%%%%%%commenting out, 12/10/09%%%%%%%%%%%%%%%%%%%%
        %for debugging:
        cor_spot = clip(cor_spot, 50000,0);
        %temp_cor_spot = newim(size(cor_spot,1), size(cor_spot,2), size(cor_spot,3)+2);
        %temp_cor_spot(:,:,1:end-1) = cor_spot;
        %temp_local_foci_mask = bdilation(local_foci_mask);
        %temp_local_foci_mask = newim(size(local_foci_mask,1), size(local_foci_mask,2), size(local_foci_mask,3)+2);
        %temp_local_foci_mask(:,:,1:end-1) = local_foci_mask;
        
        
        if sum(and(local_foci_mask>0,local_nuc>0))    %if the local foci mask is not empty
            %apparently, dip_localminima only looks at z indices 1:end-1
            %local_spot_single = dip_localminima(-(add_z_ends(cor_spot)), add_z_ends(bdilation(local_foci_mask))>0, 3, max_int, max_size, 1);   %initial foci seed detection
            %local_spot_single = dip_localminima(-temp_cor_spot, temp_local_foci_mask>0, 3, max_int, max_size, 1);   %initial foci seed detection
            try
                local_spot_single = dip_localminima(-local_spot, local_foci_mask>0, 2, max_int, max_size, 1);
                %local_spot_single = dip_localminima(-(add_z_ends(local_spot)), add_z_ends(local_foci_mask)>0, 3, max_int, max_size, 1);   %initial foci seed detection
                %local_spot_single = local_spot_single(:,:,1:end-1);
                local_foci_mask = dip_growregions(label(local_spot_single),[],local_foci_mask, dimen, 0, 'high_first');
                % Weed out dim foci (mean Int < min_int)
                % S. Costes, September 2009
                % trying with mean Int < min_int + 0.1(max_int-min_int)
                temp_ms = measure(local_foci_mask, local_spot, 'MaxVal', [], 2);
                temp_local = msr2obj(local_foci_mask,temp_ms,'MaxVal');
                local_foci_mask = temp_local>min_int;
            catch
                local_spot_single = newim(size(local_foci_mask));
            end
        else
            local_spot_single = newim(size(local_foci_mask));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %make it a binary        
        local_foci_mask = local_foci_mask>0;

        
        local_spot_single = (local_spot_single*local_nuc*local_foci_mask)>0;            %isolate seeds only within the nuclear segmented area/volume, also in the local_foci_mask
        local_label_spot_single = label(local_spot_single,2);
        %need to deal with grouped seeds
        local_label_spot_single_ms = measure(local_label_spot_single, [], {'size', 'center'}, [], 2);  %take measurement of size and center of labeled local_spot_single
        if isempty(local_label_spot_single_ms)
            large_seeds = [];
        else
            large_seeds = find(local_label_spot_single_ms.size>1);      %find spots > 1 in size
        end
        if ~isempty(large_seeds)        %if the large seed list isnt' empty
            for m = large_seeds         %for each id associated with a seed > 1 in size
                new_large_seed_coords = round(local_label_spot_single_ms(m).center);    %grab the center coordinates and round it
                local_label_spot_single(find(local_label_spot_single==m)) = 0;          %set the original pixels to 0
                local_label_spot_single(new_large_seed_coords(1), new_large_seed_coords(2)) = m;      %set the new rounded coordinate to the id
            end
        end


        local_spot_single = local_label_spot_single>0;      %the labeled version has been filtered for seeds > 1. take the binary and set it back to local_spot_single


        %grow regions back
        %local_label_foci_mask = dip_growregions(local_label_spot_single, [], local_foci_mask, 3, 0, 'high_first', 'default');
        local_label_foci_mask = dip_growregions(local_label_spot_single, [], local_foci_mask, 2, 0, 'high_first');

        %OK TO HERE%
        
        %%%added 12/20/07 - test gravity method of obtaining spot_mask.
        %%for each of the foci (in 3D), find the center slice (based off of
        %%ms.Center) and get rid of all the other slices for this foci.
        %%this "trimmed" down version of the local labeled foci mask is
        %%what is used.
        trim_local_label_foci_mask = local_label_foci_mask;
        
       
        
        %HERE WE HAVE THE LOCAL LABEL FOCI MASK, THE TRIMMED LOCAL LABEL
        %FOCI MASK AND THE SPOT SINGLE MASK. CHECK WHICH FOCI POSITIONS ARE
        %BETTER BY OVERLAYING WITH THE LOCAL SPOT IMAGE
        %%%



        %grab foci measurements
        if (max(local_label_foci_mask) == 0)
            count_spot(nuc_num) = 0;
            info_spot = [info_spot; nuc_num, NaN, NaN, NaN, NaN];
        else
            %take measurements for foci information
            foci_ms = measure(local_label_foci_mask,local_spot,{'Mean','MaxVal','Size'}, [], 2);
            spot_mean = foci_ms.mean;
            spot_max = foci_ms.maxval;
            spot_size = foci_ms.size;
            spot_ID = foci_ms.ID;
            num_spot = length(spot_mean);
            info_spot = [info_spot; ones(num_spot,1)*nuc_num, spot_ID', spot_size', spot_mean', spot_max'];
            count_spot(nuc_num) = num_spot;
        end

        %patch local_label_foci_mask into foci_mask
        foci_mask(coords(1,1):coords(1,2),coords(2,1):coords(2,2)) = foci_mask(coords(1,1):coords(1,2),coords(2,1):coords(2,2)) + local_label_foci_mask;

        %patch trim_local_label_foci_mask into trim_foci_mask
        trim_foci_mask(coords(1,1):coords(1,2),coords(2,1):coords(2,2)) = trim_foci_mask(coords(1,1):coords(1,2),coords(2,1):coords(2,2)) + trim_local_label_foci_mask;

        %patch local_spot_single into spot_single
        spot_single(coords(1,1):coords(1,2),coords(2,1):coords(2,2)) = spot_single(coords(1,1):coords(1,2),coords(2,1):coords(2,2)) + local_spot_single;

    end     %end if dimen == 3
%    return_wps_arr{nuc_num} = wps;
return_wps_arr{nuc_num} = NaN;
end     %end for nuc_num = 1:total_nuc

%fprintf('\n');
%close(h);

spot_mask = spot_single;
full_spot = foci_mask;
trim_full_spot = trim_foci_mask;
if ~exist('nuc_coords')
    nuc_coords = [];
end
if ~exist('count_spot')
    count_spot = 0;
end
spot_struct = struct(...
    'nucID',(1:total_nuc)',...
    'nucID_info','ID for each analyzed nucleus',...
    'nuc_coords',nuc_coords,...
    'nuc_coords_info','coordinates encompassing nuclei:coords(nucI,1,1):coords(nucI,1,2),coords(nucI,2,1):coords(nucI,2,2),coords(nucI,3,1):coords(nucI,3,2))',...
    'count',count_spot',...
    'count_info','number of foci/nucleus',...
    'foci',info_spot,...
    'foci_info','[nucID, fociID, foci_size, foci_meanI, foci_maxI]');
end

function new_img = add_z_ends(img)
%will add blank z "bookends" onto a stack
[w,h,d] = size(img);
new_img = newim([w, h, d+2]);
new_img(:,:,1:end-1) = img;
end