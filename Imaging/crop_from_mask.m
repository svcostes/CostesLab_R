function [img,coords] = crop_from_mask(img1,mask, boundary);
% [img,coords] = crop_from_mask(img1,mask,boundary)
% INPUT
% img is the crop of img1 delimited by mask
% boundary is the additional edge around the strict crop to be used when
% creating the final crop; Default is 1;
%
% OUTPUT
% coords are the coordinates defining the cropped region
% img = img1(coords(1,1):coords(1,2),coords(2,1):coords(2,2)); in 2D
% img = img1(coords(1,1):coords(1,2),coords(2,1):coords(2,2),coords(3,1):coords(3,2)); in 3D
% 
% See also: insert_crop, reverse of this function
%
% Sylvain Costes, July 2005, LBNL

if ~exist('boundary')
    boundary = 1;
end

if ~isempty('boundary')
    boundary = 1;
end

if min(size(img1{1}) == size(mask)) == 0 % The {1} will not crash on a color image
    error('gray image and binary image must have same dimension');
end

dimen = length(size(img1{1}));

if dimen<3
    x = sum(mask,[],2)>0;
    y = sum(mask,[],1)>0;
    left_x = max(min(find(x>0))-boundary,0);
%     if (left_x <0)
%         left_x = 0;
%     end
%     if (left_x >= size(img1,1))
%         left_x = size(img1,1)-1;
%     end
    
    right_x = min(max(find(x>0))+boundary,size(img1{1},1)-1);
%     if (right_x <0)
%         right_x = 0;
%     end
%     if (right_x >= size(img1,1))
%         right_x = size(img1,1)-1;
%     end
    
    top_y = max(min(find(y>0))-boundary,0);
%     if (top_y <0)
%         top_y = 0;
%     end
%     if (top_y >= size(img1,2))
%         top_y = size(img1,2)-1;
%     end
    
    bot_y = min(max(find(y>0))+boundary,size(img1{1},2)-1);
%     if (bot_y <0)
%         bot_y = 0;
%     end
%     if (bot_y >= size(img1,2))
%         bot_y = size(img1,2)-1;
%     end
    
    img = img1(left_x:right_x,top_y:bot_y);
    coords = [left_x,right_x;top_y,bot_y];
else
    z = sum(mask,[],3)>0;
    x = sum(z,[],2)>0;
    y = sum(z,[],1)>0;
    z = sum(sum(mask,[],1),[],2)>0;
    
    left_x = max(min(find(x>0))-boundary,0);
    right_x = min(max(find(x>0))+boundary,size(img1{1},1)-1);
    top_y = max(min(find(y>0))-boundary,0);
    bot_y = min(max(find(y>0))+boundary,size(img1{1},2)-1);
    min_z = max(min(find(z>0))-boundary,0);
    max_z = min(max(find(z>0))+boundary,size(img1{1},3)-1);    
%     left_x = min(find(x>0));
%     right_x = max(find(x>0));
%     top_y = min(find(y>0));
%     bot_y = max(find(y>0));
%     min_z = min(find(z>0));
%     max_z = max(find(z>0));
    img = img1(left_x:right_x,top_y:bot_y,min_z:max_z);
    coords = [left_x,right_x;top_y,bot_y;min_z,max_z];
end