% October 2015 - compute averages per well - If multiple processed file,
% concatenate the latest directory
% concatenate_well(nplate,well_dir)
% if nplate is the name of the directory, assumes it has the format
% P200_123433434, will store in database P200 and save summary file in
% directory.
% if nplate is a number, this assumes the plate is sprintf('P%d',nplate)
% S. Costes, Exogen Bio, 2016
% Modified on 12/10/2018 for working off mongoDB. Need to update with a way
% to tell which well is which... Later.
% S. Costes

function concatenate_well_Ames(plate_name, plate_dir)

% averaged data
% num_nuc avg_nfoci std_nfoci num_nuc_no_outl
% avg_foci_no_outl std_foci_no_outl avg_fitc avg_fitc_bg avg_nuc_ara
% avg_p2a
% averaged strings
% Well # Process_dir
fop_outf = fullfile(plate_dir, ['average_well_' plate_name '.txt'])
fop = fopen(fop_outf,'w'); % Average file for all wells
fcop = fopen(fullfile(plate_dir, ['concat_plate_' plate_name '.txt']),'w');% Concatenated output

fprintf(fop,'Well\tDirectory\tnum_nuc\tavg_nfoci\tstd_nfoci\tnum_nuc_no_outl\tavg_foci_no_outl\tstd_foci_no_outl\tavg_fitc\tavg_fitc_bg\tavg_nuc_area\tavg_nuc_dapi\tavg_p2a\tspot_fitc_sum\n');
fprintf(fcop,'Image\tDirectory\tWell#\tfitc_background\tnuc_dapi\tnuc_area\tp2a\tnfoci\tfitc_mean\tfitc_sum\n');

% Read process directory for each well. If multiple process directory, read
% latest analysis directory only
list_well = dir(plate_dir);
list_well = list_well(~ismember({list_well.name},{'.','..','Color_images'})) %ignore . and .. files as well as Color_images folder.
list_well = list_well([list_well.isdir]) % keep only
num_well = length(list_well);
for i_well = 1:num_well
    disp(i_well)
    well_dirs = dir(fullfile(plate_dir, list_well(i_well).name, 'Foci_Processed*'));
    disp(fullfile(plate_dir, list_well(i_well).name, 'Foci_Processed*'))
    num_process = length(well_dirs);
    if num_process == 0
        continue
    else if num_process == 1
            well_file = dir(fullfile(plate_dir, list_well(i_well).name,  well_dirs(1).name, 'nuc_summary_foci_analysis*'));
            keep_ind = 1;
        else if num_process >1 % well file to most recent based on file datenum
                [tmp,keep_ind] = sort([well_dirs.datenum]);
                well_file = dir(fullfile(plate_dir, list_well(i_well).name, well_dirs(keep_ind(end)).name, 'nuc_summary_foci_analysis*'));
            end
        end
        disp(well_dirs)
        disp(list_well(i_well).name)
        disp(well_file)
        well_file = well_file(1).name;
        % find corresponding treatment for each well
        well_str = list_well(i_well).name;
        filename = fullfile(plate_dir, list_well(i_well).name, well_dirs(keep_ind(end)).name, well_file);
        fprintf('Reading %s\n',filename);
        % Index (data or text as appropriate)
        % 1 nucID
        % 2 image
        % 3 fitc_background
        % 4 txred_background
        % 5 nuc_dapi
        % 6 nuc_area
        % 7 p2a
        % 8 nfoci
        % 9 spot_mean (NAN if no foci - use string)
        % 10 fitc_mean
        % 11 fitc_sum
        % 12 txred_mean
        % 13 txred_sum
        % 14 foci_dapi - skip this one.
        fid = fopen(filename);
        datat = textscan(fid,'%f%s%f%f%f%f%f%f%s%f%f%f%f%*[^\n]','headerlines',1,'delimiter','\t');
        fclose(fid);
        num_nuc = length(datat{1});
        for i_nuc = 1:num_nuc
            fprintf(fcop,'%s\t%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',datat{2}{i_nuc},list_well(i_well).name,well_dirs(keep_ind(end)).name,...
                datat{3}(i_nuc),datat{5}(i_nuc),datat{6}(i_nuc),datat{7}(i_nuc),datat{8}(i_nuc),datat{10}(i_nuc),datat{11}(i_nuc));
        end
        p2a = datat{7};
        nfoci = datat{8};
        spot_mean = datat{9};
        nan_ind = find(strcmp(spot_mean,'NaN'));
        spot_mean = str2double(datat{9});
        spot_mean(nan_ind) = 0;
        nuca = datat{6};
        fitc = datat{10};
        % Find index of data to keep within standard deviation
        ki = find(and(and(p2a<(mean(p2a)+2*std(p2a)),and(nuca<mean(nuca)+1.5*std(nuca),nuca>mean(nuca)-1.5*std(nuca))),fitc<(mean(fitc)+2*std(fitc))));
        kept_nuc = length(ki);
        fprintf(fop,'%s\t%s\t%d\t%5.2f\t%5.2f\t%d\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\n',list_well(i_well).name,well_dirs(keep_ind(end)).name,...
            num_nuc,mean(datat{8}),std(datat{8}),kept_nuc,mean(nfoci(ki)),std(nfoci(ki)),mean(fitc(ki)),mean(datat{3}(ki)),mean(nuca(ki)),mean(datat{5}(ki)),mean(p2a(ki)),mean(nfoci(ki).*spot_mean(ki)));
        %catch ME
        %    fprintf('Problem reading %s\n', [well_dir list_well(i_well).name]);
    end
end
fclose(fop);
fclose(fcop);
end
