function [img] = insert_crop(crop_img,img,coords,mode);
% [img] = insert_crop(crop_img,img,coords,mode)
% put crop_img into img at coordinates coords
% mode: 'add', will add to what is in img (default). 'new', will overwrite
% anything where images is insterted.
%
% Format for coords
% [coordsX;
%  coordsY;
%  coordsZ]
% coordsX is either one or two scalar vector. If one, assume upper left
% coordinates are passed. If two, upper left and lower right coordinates
% are passed.
%
% crop_img = img(coords(1,1):coords(1,2),coords(2,1):coords(2,2)); in 2D
% crop_img = img(coords(1,1):coords(1,2),coords(2,1):coords(2,2),coords(3,1):coords(3,2)); in 3D
%
% See also: crop_from_mask, reverse of this function
%
% Sylvain Costes, October 2009, LBNL

if sum(size(crop_img)>size(img))
    error('Cannot put a larger image inside a smaller image');
end

if ~exist('mode','var')
    mode = 'add';
end

switch(lower(mode))
    case 'add'
        if size(coords,2) > 1
            if length(size(img)) > 2
                img(coords(1,1):coords(1,2),coords(2,1):coords(2,2),coords(3,1):coords(3,2)) = crop_img + img(coords(1,1):coords(1,2),coords(2,1):coords(2,2),coords(3,1):coords(3,2));
            else
                img(coords(1,1):coords(1,2),coords(2,1):coords(2,2)) = crop_img + img(coords(1,1):coords(1,2),coords(2,1):coords(2,2));
            end
        else
            xmax = coords(1,1) + size(crop_img,1)-1;
            ymax = coords(2,1) + size(crop_img,2)-1;
            if length(size(img)) > 2
                zmax = coords(3,1) + size(crop_img,3)-1;
                img(coords(1,1):xmax,coords(2,1):ymax,coords(3,1):zmax) = crop_img + img(coords(1,1):xmax,coords(2,1):ymax,coords(3,1):zmax);
            else
                img(coords(1,1):xmax,coords(2,1):ymax) = crop_img + img(coords(1,1):xmax,coords(2,1):ymax);
            end
        end
    case 'new'
        if size(coords,2) > 1
            if length(size(img)) > 2
                img(coords(1,1):coords(1,2),coords(2,1):coords(2,2),coords(3,1):coords(3,2)) = crop_img ;
            else
                img(coords(1,1):coords(1,2),coords(2,1):coords(2,2)) = crop_img;
            end
        else
            xmax = coords(1,1) + size(crop_img,1)-1;
            ymax = coords(2,1) + size(crop_img,2)-1;
            if length(size(img)) > 2
                zmax = coords(3,1) + size(crop_img,3)-1;
                img(coords(1,1):xmax,coords(2,1):ymax,coords(3,1):zmax) = crop_img;
            else
                img(coords(1,1):xmax,coords(2,1):ymax) = crop_img;
            end
        end
end