function [final_arr, wavelet_planes] = james_a_trous(img, level)
%This function will compute the a trous wavelet transformation on img, up
%to the number of levels specified in level. It will return each successive
%transformation in a stacked array final_arr, which includes the original
%image as the first slice. wavelet_planes returns the difference between
%the transformations (final_arr has n slices, wavelet_planes will have n-1
%planes).
%
%To recontruct the original image, sum up the wavelet_planes in addition to
%the last transformation slice i.e. reconstructed_img = final_arr(:,:,n) +
%wavelet_planes(:,:,1) + wavelet_planes(:,:,2) + ... +
%wavelet_planes(:,:,n-1);
%
%ex: img = readim('test.tif');
%    [final_arr, wps] = james_a_trous(img, 4);
%
%NOTE: img must be stretched to 8 bit before passing in. (ex:
%   (img-250)/(3000-250)*255;)
%
%James Chen, jlchen@lbl.gov, 07/12/07

    dimen = length(size(img));      %determine if 2d or 3d    
    
    
    img_arr = double(img);      %convert to matrix of numbers
    [w, h, d] = size(img);      %get dimensions
    
    if dimen == 2
        final_arr(:,:,1) = img_arr;     %load first final_arr with original image
    else
        final_arr(:,:,:,1) = img_arr;   %load first final_arr with original stack
    end
    
    kernal_base = [1/16 1/4 3/8 1/4 1/16];

    for i = 1:level     %for each level, starting at 1
%         tic;
%         fprintf('LEVEL: %d\n', i);
        row_arr = [];       %this variable is overwritten with different dimensions each iteration. must clear beforehand
        col_arr = [];       %this variable is overwritten with different dimensions each iteration. must clear beforehand
        num_of_zeroes = 2^(i-1)-1;      %number of zeroes inbetween kernal_base elements. depends on level.
        %num_of_zeroes = i-1;
        kernal_length = 4*num_of_zeroes + 5;    %total length of kernal, including the above zeroes
        kernal = zeros(1,kernal_length);    %kernal is initially zeroes
        kernal_counter = 1;         %kernal counter for indexing through kernal_base

        for j = 1:(num_of_zeroes+1):kernal_length
            kernal(j) = kernal_base(kernal_counter);        %evenly divide the kernal base elements into kernal
            kernal_counter = kernal_counter +1;
        end        

        %kernal vector ex:
        %level 1: [1/16 1/4 3/8 1/4 1/16]
        %level 2: [1/16 0 1/4 0 3/8 0 1/4 0 1/16];
        %level 3: [1/16 0 0 0 1/4 0 0 0 3/8 0 0 0 1/4 0 0 0 1/16];
        %the number of zeroes at each level follows 2^(i-1)-1
        
        edge_buffer = kernal_length;
        
        %we now have our kernal vector, for level i
        dead_zone = (kernal_length + 1)/2;      %amount we need to clip from the beginning and end of the convolved image. extra length resulting from convolution;

        if dimen == 2
            
            %create an extended version of img_arr, which replicates the
            %end bounding rows and columns as a buffer zone.
            new_img_arr = create_extended_img_arr(img_arr, edge_buffer);
                       
            %2d image
            %1st approach: separable convolutions
            temp = separable_2d_conv(new_img_arr, kernal, dead_zone);       %perform convolution on extended array
            final_arr(:,:,i+1) = temp(edge_buffer+1:end-edge_buffer, edge_buffer+1:end-edge_buffer);        %trim array down to original size, store in final_arr

            wavelet_planes(:,:,i) = final_arr(:,:,i)-final_arr(:,:,i+1);            %compile a stack of the wavelet planes (the difference between the transformed images, beteween scales)
            
            %%2nd approach: conv2
            %
            %final_arr(:,:,i+1) = conv2(img_arr, kron(kernal,kernal'),'same');
            %wavelet_planes(:,:,i) = final_arr(:,:,i) - final_arr(:,:,i+1);
            %%end 2nd approach
        else

            %1st approach: separable convolutions
            for u = 1:size(img,3)
                new_img_arr(:,:,u) = create_extended_img_arr(img_arr(:,:,u), edge_buffer);
            end

     
            temp = separable_3d_conv(new_img_arr, kernal, dead_zone);   %compile a stack of the transformed image at each scale
    
       
            final_arr(:,:,:,i+1) = temp(edge_buffer+1:end-edge_buffer, edge_buffer+1:end-edge_buffer,:);            
            wavelet_planes(:,:,:,i) = final_arr(:,:,:,i)-final_arr(:,:,:,i+1);      %compile a stack of the wavelet planes (the difference between the transformed images, beteween scales)
        
            %%2nd approach: convn
            %kernal = kron3d(kernal);
            %final_arr(:,:,:,i+1) = convn(img, kernal, 'same');
            %wavelet_planes(:,:,:,i) = final_arr(:,:,:,i)-final_arr(:,:,:,i+1);
            %%end 2nd approach
        end
        clear new_img_arr;
 %       t1 = toc;
%         fprintf('level time: %3.3f\n', t1);
    end     %end for i = 1:level

    
end     %end function james_a_trous

%------------------------------------------------------------------------------------------------------------------------
function new_img_arr = create_extended_img_arr(img_arr, edge_buffer)
    %extends img_arr which is x by y, to x+length by y+length. Mirrors the
    %image to edge_buffer extra region.
    
    %create new_img_arr with additional boundary regions
    new_img_arr = double(newim(size(img_arr,2) + 2*edge_buffer, size(img_arr,1) + 2*edge_buffer));

    
    
    %copy img_arr into center of new_img_arr
    new_img_arr(edge_buffer+1:end-edge_buffer, edge_buffer+1:end-edge_buffer) = img_arr;
    

    [h, w] = size(img_arr);
    
    if edge_buffer >= w
        r_start = edge_buffer-w + 1;
    else
        r_start = 1;
    end
    
     
    %fill left boundary region
    for r = r_start:edge_buffer
        new_img_arr(edge_buffer+1:end-edge_buffer,r) = img_arr(:,edge_buffer-r+1);
    end

    %fill right boundary region
    for r = r_start:edge_buffer
        new_img_arr(edge_buffer+1:end-edge_buffer,end-r+1) = img_arr(:,end-edge_buffer+r);
    end

    if edge_buffer >= h
        r_start = edge_buffer-h + 1;
    else
        r_start = 1;
    end
    
    %fill top boundary region
    for r = r_start:edge_buffer
        new_img_arr(r,1:end) = new_img_arr(2*edge_buffer-r+1,1:end);
    end

    %fill bottom boundary region
    for r = r_start:edge_buffer
        new_img_arr(end-r+1,1:end) = new_img_arr(end-2*edge_buffer+r, 1:end);
    end

end     %end function new_img_arr


%------------------------------------------------------------------------------------------------------------------------

function z_arr = separable_3d_conv(img_arr, kernal, dead_zone, z_option)
    %performs a 3d convolution as a combination of 1d combinations along
    %each dimension. if z_option is 1, do z-axis convolution.
    if ~exist('z_option')
        z_option = 0;
    end
    [h, w, d] = size(img_arr);      %size(img_arr) and size(img) are transposed
    for j = 1:d     %convolve each slice as a 2d image
        temp_arr(:,:,j) = separable_2d_conv(img_arr(:,:,j), kernal, dead_zone);
    end     %end for j = 1:d

    if (z_option)
        for j = 1:h   %convolve each z-stack column with the kernal
            for k = 1:w
                z_arr(j,k,:) = conv(kernal, reshape(temp_arr(j,k,:),1,d,1));       %temp_arr is transposed from temp. thus we must access by (h,w) instead of (w,h)
            end     %end for k = 1:h
        end     %end for j = 1:w

        z_arr = z_arr(:,:,dead_zone:end-dead_zone+1);        %clip back to original size of image.
    else
        z_arr = temp_arr;
    end
    %z_arr has the compiled 3d convolution array

end

%------------------------------------------------------------------------------------------------------------------------

function col_arr = separable_2d_conv(img_arr, kernal, dead_zone)
    %performs a 2d convolution as a combination of 1d combinations along
    %each dimension.

    [h, w] = size(img_arr);     %size(img_arr) and size(img) are transposed
    %do row by row convolution of the image with the kernal
    for j = 1:h
        row_arr(j,:) = conv(kernal, img_arr(j,:));  %row i
    end
    %row_arr now has dimensions [w+kernal_length-1, h]
    
    row_arr = row_arr(:,dead_zone:end-dead_zone+1);     %clip back to original size of image.
    
    %do col by col convolution of the previously convolved image, with the kernal
    for j = 1:w
        col_arr(:,j) = conv(kernal,row_arr(:,j));       %col i
    end
    
    col_arr = col_arr(dead_zone:end-dead_zone+1,:);   %clip back to original size of image.
    %col_arr has the compiled 2d convolution array

end     %end separable_2d_conv

%------------------------------------------------------------------------------------------------------------------------

function k3d = kron3d(kernal)
    k2d = kron(kernal,kernal');
    for i = 1:length(kernal)
        k3d(:,:,i) = k2d*kernal(i);
    end
end     %end function kron3d(kernal)