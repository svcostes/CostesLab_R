function struct_print(ofp, thestruct, header_opt, PRINT_FLAG)
% struct_print(ofp, thestruct, header_opt)
% Will print content of thestruct to output file ofp. If header_opt>1 then
% header are printed first.
% Print done horizontally
%
% This function assumes that each vectors for each field has the same
% length and are properties of the same object. Thus each row corresponds
% to the properties reported in the structure variable and listed as a
% given field.
%
% Charlie Chen, LBNL, Feb, 2011
% edited for supporting printing string in the struct
%
% Sylvain Costes, LBNL, October 2009

if nargin < 3
    header_opt = 0;
end

if ~exist('PRINT_FLAG') || isempty(PRINT_FLAG); PRINT_FLAG = false; end

try
    fn = fieldnames(thestruct);

    if header_opt % Print header if option is on
        for n=1:length(fn)
            fprintf(ofp,'%s\t',fn{n});
        end
        fprintf(ofp,'\n');
    end

    c = struct2cell(thestruct);

    find_size_flag = 1;
    cnt = 1;
    num_obj = 1; %default, in case cannot find numerical arrays
    while find_size_flag
        if ischar(c{cnt})
            cnt = cnt + 1;
        else
            num_obj = length(c{cnt});
            find_size_flag = 0;
        end
    end

    for nf = 1:num_obj % Loop on # of objects
        for n = 1:length(fn) % Loop on all property values for one object (one row output)
            try
                if ischar(c{n}) % if input is a string
                    fprintf(ofp,'%s\t',c{n});
                elseif isinteger(c{n}(nf))
                    fprintf(ofp,'%d\t',c{n}(nf));
                elseif isnumeric(c{n}(nf))
                    fprintf(ofp,'%f\t',c{n}(nf));
                    %fprintf('%f\t',c{n}(nf));
                else
                    fprintf(ofp,'%s\t',c{n}{nf});
                    %fprintf('%s\t',c{n}{nf});
                end
            catch
                fprintf(ofp,'NaN\t');
                %fprintf('NaN\t');
            end
        end
        fprintf(ofp,'\n');
        %fprintf('\n');
    end
catch
    if PRINT_FLAG
        fprintf('Could not print anything for this measurement file\n');
    end
end
