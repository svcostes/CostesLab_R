function img_out = map(f,img_3D,params)
% Map function useful to repeat in 3D a 2D operation
% img_out = map(f,img_3D,params)
%
% Example: map('tophat',img,'7');
% Will apply tophat for each 2D slice such that img_out(:,:,i_slice) =
% tophat(img(:,:,i_slice),7);
%
% Sylvain Costes, LBNL, Oct. 2010
%
num_slice = size(img_3D,3);
img_out = newim(size(img_3D));
for k = 1:num_slice
    if exist('params')
        img_out(:,:,k-1) = eval([f '(squeeze(img_3D(:,:,k-1)),' params ')']);
    else
        img_out(:,:,k-1) = eval([f '(squeeze(img_3D(:,:,k-1)))']);
    end
end