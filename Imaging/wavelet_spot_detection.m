function [foci_mask, return_wps] = wavelet_spot_detection(img, mask, levels, k)
%attempt to do spot detection using james_a_atrous wavelet transformation
%for 2d image.
%
%ex: foci_mask = wavelet_spot_detection(readim('test.tif'), 3,3);
%ex: foci_mask = wavelet_spot_detection(readim('test3d.ics'),3,4); 
%

    dimen = length(size(img));              %find dimensions of image
    [final_arr, wps] = james_a_trous(img,levels);       %do a trous transformations

    [w,h,d] = size(img);      %grab dimensions of the original image
    ld = 1;                 %threshold for correlation matrix
    %k = 3;                  %constant for wavelet plane threshold
   
    for i = 1:levels        %for each level
%       fprintf('wave LEVEL: %d\n', i);
       temp_1d = [];       %clear temp_1d
       temp_image = [];       %clear temp_image

       if dimen == 2
           temp_image = wps(:,:,i);   %copy current wavelet_plane into temp_2d
       else
           temp_image = wps(:,:,:,i);
       end
        %temp_1d = reshape(temp_image, 1, w*h*d);     %reshape into 1d
        
       
        temp_1d = temp_image(:);


        
        sigma_i = (medmad(temp_1d))/0.67;       %find the median deviation
        wp_thresh = k*sigma_i;                  %wavelet plane threshold for this level

        temp_image(find(temp_image < wp_thresh)) = 0;     %set pixels less than this threshold in the wavelet plane to 0

      
        
       if dimen == 2
           wps2(:,:,i) = temp_image;                  %transfer it to a stacked 3d variable wps2
       else
           wps2(:,:,:,i) = temp_image;
       end

    end     %end for i = 1:levels
   

    if dimen == 2
        cor_mat = wps2(:,:,1).*wps2(:,:,2);         %multiply the first two wavelet planes, element by element
        for i = 3:levels                            %do the product with the remaining wavelet planes
            cor_mat = cor_mat.*wps2(:,:,i);
        end     %end for i = 3:levels
    else
        cor_mat = wps2(:,:,:,1).*wps2(:,:,:,2);
        for i = 3:levels
            cor_mat = cor_mat.*wps2(:,:,:,i);
        end
    end

    abs_cor_mat = abs(cor_mat);                 %take the absolute value of the final correlation matrix product
    cor_mat(find(abs_cor_mat < ld)) = 0;        %set pixels less than correlation matrix threshold to 0

 %   foci_mask = fill_in_holes((dip_image(cor_mat>0)*mask)>0);          %this is the final spot detected image
     foci_mask = label((dip_image(cor_mat>0)*bdilation(mask))>0,dimen,3,100000)>0;          %this is the final spot detected image
   % Added dilation 9/10/10, SVC, so that border foci are not rejected due
   % to bad nuc seg. Can eliminate foci after, not here...
    
%     trim_foci_mask = newim(size(foci_mask));
%     trim_foci_spot = newim(size(foci_mask));
%     %find optimal z slice for each foci based on measurement. return as trim_foci_mask
%     if max(foci_mask ~= 0)  %if foci_mask is empty
%         foci_mask_ms = measure(dip_image(uint8(foci_mask)), [], 'center', [], 2);     %get ids of all foci, and centers
%         for k = 0:size(foci_mask,3)-1      %for each slice
%             fprintf('k: %d\n', k);
%             indices = find(round(foci_mask_ms.Center(3,:))==k);    %
%             keep_ids = foci_mask_ms.ID(indices);
%             if ~isempty(keep_ids)
% 
%                 %HERE NEED TO FIND A WAY TO WEED OUT BAD
%                 %SLICES THAT DON'T WANT.
%                 temp_cff = double(foci_mask(:,:,k));
%                 keep_coords = compare2Arrays(temp_cff(:), keep_ids);
%                 temp_trim_foci_mask = newim(size(squeeze(foci_mask(:,:,k))));
%                 temp_trim_foci_mask(keep_coords-1) = 1;
%                 temp_trim_foci_mask = temp_trim_foci_mask*squeeze(foci_mask(:,:,k));
%                 
%                 %foci_spot_center = measure(dip_image(uint8(temp_trim_foci_mask)), squeeze(img(:,:,k)), 'gravity', []);
%                 
%                 trim_foci_mask(:,:,k) = temp_trim_foci_mask;
%             end
%         end
%     else
%         trim_foci_mask = foci_mask;
%     end
    
    
    
%     crop_fitc_foci_ms = measure(dip_image(uint8(crop_fitc_foci)), [], 'center', [], 2);     %get ids of all foci, and centers
%     trim_crop_fitc_foci = newim(size(crop_fitc_foci));
%     trim_crop_fitc_spot = newim(size(crop_fitc_foci));
%     [h w] = size(trim_crop_fitc_spot);
%     for k = 0:size(crop_fitc_foci,3)-1      %for each slice
%         fprintf('k: %d\n', k);
%         indices = find(round(crop_fitc_foci_ms.Center(3,:))==k);    %
%         keep_ids = crop_fitc_foci_ms.ID(indices);
%         if ~isempty(keep_ids)
% 
%             %HERE NEED TO FIND A WAY TO WEED OUT BAD
%             %SLICES THAT DON'T WANT.
%             temp_cff = double(crop_fitc_foci(:,:,k));
%             keep_coords = compare2Arrays(temp_cff(:), keep_ids);
%             temp_trim_crop_fitc_foci = newim(size(squeeze(crop_fitc_foci(:,:,k))));
%             temp_trim_crop_fitc_foci(keep_coords-1) = 1;
%             temp_trim_crop_fitc_foci = temp_trim_crop_fitc_foci*squeeze(crop_fitc_foci(:,:,k));
% 
%             spot_center = measure(dip_image(uint8(temp_trim_crop_fitc_foci)), (squeeze(crop_fitc(:,:,k))), 'gravity', []);
%             spot_center_coords = round(spot_center.Gravity)';
%             %trim_crop_fitc_spot(
% 
%             trim_crop_fitc_spot(spot_center_coords(:,1)*w+spot_center_coords(:,2)+k*w*h) = 1;
% 
%             trim_crop_fitc_foci(:,:,k) = temp_trim_crop_fitc_foci;
%         end
%         %trim_crop_fitc_foci(:,:,k) = crop_fitc_foci(:,:,k) == keep_ids;
%     end
    
    
return_wps = wps;

end     %end function wavelet_spot_detection

%--------------------------------------------------------------------------------------------------------------

function temp_img = fill_in_holes(temp_img)
    %fill in holes. for 3d, fill in holes for each SLICE independently.
    if (length(size(temp_img)) == 2)
        temp_img = label(temp_img, 2, 2)>0;             %get rid of single pixel objects in the slice.
        temp_img = berosion(bdilation(temp_img,1,2),1,2);       %erode and dilate. to join sparsely separated objects
        ms = measure(label(~temp_img), [], 'size');     %find size of each labelled object in the image, inverted
        bg_id = find(ms.size == max(ms.size));          %find the largest sized object in the image. should be bg
        if bg_id ~= 1                                   %if the largest id is not 1, (not the bg)
            temp = and(label(~temp_img,1) ~= 1, label(~temp_img,1) ~= bg_id);   
        else
            temp = label(~temp_img, 1) ~= 1;
        end
        temp_img = label(or(temp, temp_img))>0;
    else
        for p = 1:size(temp_img,3)
            temp_img2 = squeeze(temp_img(:,:,p-1));         %get 2d slice
            temp_img2 = label(temp_img2, 2, 2)>0;           %get rid of single pixel objects in the slice.
            temp_img2 = berosion(bdilation(temp_img2,1,2),1,2);         %erode and dilate. to join sparsely separated objects
            ms = measure(label(~temp_img2), [], 'size');     %find size of each labelled object in the image, inverted
            bg_id = find(ms.size == max(ms.size));          %find the largest sized object in the image. should be bg
            if bg_id ~= 1                                   %if the largest id is not 1, (not the bg)
                temp = and(label(~temp_img2,1) ~= 1, label(~temp_img2,1) ~= bg_id);
            else
                temp = label(~temp_img2, 1) ~= 1;
            end
            temp_img2 = label(or(temp, temp_img2))>0;
            
            temp_img(:,:,p-1) = temp_img2;              %put modified slice back into stack
        end     %end for p = i:size(temp_img,3)
    end
end     %end fill_in_holes
