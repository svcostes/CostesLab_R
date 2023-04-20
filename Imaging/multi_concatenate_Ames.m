% Concatenate all plates within a directory
% S. Costes, 12/10/2018 - NASA Ames
%
ImageDir = 'D:\';
plate_list = dir(ImageDir);
num_dir = length(plate_list);

for i_dir = 1:num_dir
    if length(plate_list(i_dir).name)>2
        concatenate_well_Ames(plate_list(i_dir).name,ImageDir);
    end
end