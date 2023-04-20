% October 2015 - compute averages per well - If multiple processed file,
% concatenate the latest directory
% concatenate_well(nplate)
% if nplate is the name of the directory, assumes it has the format
% P200_123433434, will store in database P200 and save summary file in
% directory.
% if nplate is a number, this assumes the plate is sprintf('P%d',nplate)
% S. Costes, Exogen Bio, 2016

function concatenate_well(nplate)

if ischar(nplate) % allow nplate to be a string with current dir to be analyzed
    plate_ind = strfind(nplate,'_');
    plate_name = nplate(1:plate_ind-1);
    cur_dir = ['F:\Images\' nplate '\'];
else
    plate_name=sprintf('P%d',nplate);
    cur_dir_d =dir(['F:\Images\' plate_name '*\']);% directory
    
    if length(cur_dir_d) ==1
        cur_dir= ['F:\Images\' cur_dir_d.name '\'];
    elseif length(cur_dir_d) > 1
        [tmp,keep_ind] = sort([cur_dir_d.datenum]);
        k=0;
        while k < length(cur_dir_d) & ~cur_dir_d(keep_ind(end-k)).isdir
            k=k+1;
        end
        cur_dir= ['F:\Images\' cur_dir_d(keep_ind(end-k)).name '\'];
    else
        fprintf('please launch concatenate_well manually - system can''t find the proper directory');
    end
end

import com.mongodb.*;
import org.*;
import com.mongodb.util.*;
import org.bson.types.*;

db= loadmongodb();
wc= WriteConcern(1);
well_coll= db.getCollection('wells_info');

well_sample=well_coll.findOne(BasicDBObject('plate',plate_name));
keys = well_sample.keySet.toArray;
key='';
for i_keys=2:length(keys)
    key = [key keys(i_keys) '\t'];
end
list_well = dir(cur_dir);
num_well = length(list_well);
% averaged data
% num_nuc avg_nfoci std_nfoci num_nuc_no_outl
% avg_foci_no_outl std_foci_no_outl avg_fitc avg_fitc_bg avg_nuc_ara
% avg_p2a
cat_data = [];
% averaged strings
% Well # Process_dir
cat_str = [];
fop = fopen([cur_dir 'average_well_' plate_name '.txt'],'w'); % Average file for all wells
fcop = fopen([cur_dir 'concat_plate_' plate_name '.txt'],'w');% Concatenated output

fprintf(fop,[ key 'Directory\tnum_nuc\tavg_nfoci\tstd_nfoci\tnum_nuc_no_outl\tavg_foci_no_outl\tstd_foci_no_outl\tavg_fitc\tavg_fitc_bg\tavg_nuc_area\tavg_nuc_dapi\tavg_p2a\tspot_fitc_sum\n']);
fprintf(fcop,'Image\tDirectory\tWell#\tfitc_background\tnuc_dapi\tnuc_area\tp2a\tnfoci\tfitc_mean\tfitc_sum\n');

% Read process directory for each well. If multiple process directory, read
% latest analysis directory only
for i_well = 1:num_well
    try
        well_dir = dir([cur_dir list_well(i_well).name '/Foci_Processed*']);
        num_process = length(well_dir);
        if num_process == 1
            well_file = dir([cur_dir list_well(i_well).name '/' well_dir(1).name '/nuc_summary_foci_analysis*']);
            keep_ind = 1;
        else if num_process >1
                [tmp,keep_ind] = sort([well_dir.datenum]);
                well_file = dir([cur_dir list_well(i_well).name '/' well_dir(keep_ind(end)).name '/nuc_summary_foci_analysis*']);
            end
        end
        well_file = well_file(1).name;
        % find corresponding treatment for each well
        well_str = list_well(i_well).name;
        filename = [cur_dir list_well(i_well).name '/' well_dir(keep_ind(end)).name '/' well_file];
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
            fprintf(fcop,'%s\t%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',datat{2}{i_nuc},list_well(i_well).name,well_dir(keep_ind(end)).name,...
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
        
        value='';
        for i_value = 2:length(keys)
            keyStr = char(keys(i_value));
            well_i=well_coll.findOne(BasicDBObject('name',[plate_name well_str]));
            value = [value char(well_i.get(keyStr))  '\t'];
        end
        
        
        fprintf(fop,[value '%s\t%d\t%5.2f\t%5.2f\t%d\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\n'],well_dir(keep_ind(end)).name,...
            num_nuc,mean(datat{8}),std(datat{8}),kept_nuc,mean(nfoci(ki)),std(nfoci(ki)),mean(fitc(ki)),mean(datat{3}(ki)),mean(nuca(ki)),mean(datat{5}(ki)),mean(p2a(ki)),mean(nfoci(ki).*spot_mean(ki)));
    catch ME
        fprintf('Problem reading %s\n', [cur_dir list_well(i_well).name]);
    end
end



fclose(fop);
fclose(fcop);

