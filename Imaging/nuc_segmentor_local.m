% nuc_segmentor_local, will segment image with minimum algorithm for speed.
%
% [seg_img] = nuc_segmentor_local(img, nuc_rad, ss_factor, th_level, edge_option, refine_option)
%
% input: img: image to segment
%        nuc_rad: radius of nucleus, used to eliminate smaller object and
%        help object separation. Also used to blur image to get local
%        background.
%        ss_factor: X,Y subsampling. The higher the faster but the rougher
%        contours of the nucleus
%        th_level: percentage above local background to segment. Default
%        1.1 (i.e. 10% above)
%        edge_option: 1 (default). IF positive, will removed nuclei
%        touching edges of the image
%        refine_option: 1 (default). Will refine each individual nuclei 
%        based on initial segmentation.
%
% Separation of nuc based on distance transform.
%
% October 2009, Sylvain Costes, Lawrence Berkeley Lab
% Modified on October 2012 from nuc_segmentor_lite3 to use local threshold
% instead.
%
% Example to segment a DAPI image "nuc" with average nuclear size 10, with
% subsampling the image by a factor 2 to reduce size, using default local
% threshold (1.1.), removing nuclei touching the edge and without final
% refinement (recommended - buggy):
%
% mask = nuc_segmentor_local(nuc,10,2,[],1,0);
% 


function [seg_img,max_ID] = nuc_segmentor_local(img, nuc_rad, ss_factor, th_level,edge_option,refine_option)
% Get image size
img = squeeze(img);
[w,h,d] = size(img);
%img = reshape(img,[w,h,d]); % make sure image is in 3D format for segmentation
dim = length(size(squeeze(img)));

%Subsample image only in X,Y plane
radS = round(nuc_rad/ss_factor);
if dim>2
    imgS = img(0:ss_factor:end,0:ss_factor:end,:);
else
    imgS = img(0:ss_factor:end,0:ss_factor:end);
end

% blur image with Gaussian to get local backround. Blur in 2D only
% dipsetpref('NumberOfThreads',1); % gaussf is buggy in multi thread mode
if dim>2
    blur_img = map('gaussf',imgS,num2str(radS*2));
else
    blur_img = gaussf(imgS,radS*2);
end

% Threshold based on local background
if exist('th_level','var')
    if isempty(th_level)
        maskS = imgS > 1.1*blur_img;
    elseif isnan(th_level)
        maskS = imgS > 1.1*blur_img;
    else
        maskS = imgS > th_level*blur_img;
    end
else
     maskS = imgS > 1.1*blur_img;
end

% Set edge_option to 1 if not entered. Default
if ~exist('edge_option','var')
    edge_option =1;
end

if ~exist('refine_option','var')
    refine_option =1;
end

if edge_option
    if dim>2 % Get rid of any object on edges
        maskS(:,0,:) = 1; % Will get rid of anything on edge
        maskS(0,:,:) = 1;
        maskS(:,end,:) = 1;
        maskS(end,:,:) = 1;
    else
        maskS(:,0) = 1; % Will get rid of anything on edge
        maskS(0,:) = 1;
        maskS(:,end) = 1;
        maskS(end,:) = 1;
    end
    maskS = label(maskS>0) > 1;
end

% This will fill holes in mask
if dim>2
    for i=0:d-1
        %    maskS(:,:,i) = berosion(bdilation(maskS(:,:,i)>0,2),3); % To make a tighter fit around the nucleus. May lose some foci though.
        inv_mask = label(~maskS(:,:,i));
        ms = measure(inv_mask,[],'size');
        [max_size,max_ind] = max(ms.size);
        max_ID = ms.ID(max_ind);
        maskS(:,:,i) = ~(inv_mask == max_ID); % Fill holes. Done for each slice, it works better.
    end
else
        inv_mask = label(~maskS);
        ms = measure(inv_mask,[],'size');
        [max_size,max_ind] = max(ms.size);
        max_ID = ms.ID(max_ind);
        maskS = ~(inv_mask == max_ID); % Fill holes. Done for each slice, it works better.
end
maskS = label(maskS,dim,round((pi*radS^dim)/2),1e6); % Remove small objects to clean up before spliting remaining
num_nuc = max(maskS);
maskR = newim(size(maskS)); 
if refine_option>0 % In case of flag on, refine each individal nuclei
    fprintf('Refining segmentation on individual nuclei\n');
    for i_nuc = 1 :num_nuc
        [crop_nuc,Ccoords] = crop_from_mask(imgS,maskS==i_nuc,10);
        crop_mask = label(threshold(crop_nuc));
        ms = measure(crop_mask,[],'size'); % Need to keep the biggest object as the nucleus
        [max_size,max_ind] = max(ms.size);
        crop_mask = crop_mask == ms.ID(max_ind);
        maskR = insert_crop(crop_mask,maskR,Ccoords,'add');
    end
else
    maskR = maskS;
end

% Use distance transform to separate touching nuclei
%fprintf('Split object\n');
maskS = object_separation_2D(maskR>0,radS,imgS);
seg_img = maskS;
% % Extrapolate back nuclear image
% if ss_factor>1
%     fprintf('Interpolating mask image back to full size\n');
%     if dim>2
%         for j = 1:d
%             seg_img = squeeze(oversample(maskS,[w,h,d]));
%             % seg_img(:,:,j-1) = resample(seg_img(:,:,j-1),1/ss_factor);
%         end
%     else
%         seg_img = squeeze(oversample(maskS,[w,h]));
%         % seg_img = resample(maskS,1/ss_factor);
%     end
% else
%     seg_img = maskS;
% end
% seg_img = dip_image(uint16(seg_img));