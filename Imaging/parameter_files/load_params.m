function params = load_params(param_file)
% load_params(param_file) will load input information for processing images
% in a directory for various spot detection processes.
% If given blank file or inaccurate path, it will prompt for this
% information in a dialog box
%
% Sylvain Costes and James Chen, September 2009, LBNL
%

if isunix %if unix, make sure all slash are in the right direction or it wont find the file...
    slash_str = '/';
else
    slash_str = '\';
end

try
    try
        param_import = textread(param_file,'%s');
    catch
        [param_file,param_dir] = uigetfile('*.txt','Choose test_params.txt file'); % Option to locate file. If cancel, it will crash and first try will be active
        param_import = textread([param_dir slash_str param_file], '%s');
    end
    params.fitc_index = [];
    params.txred_index = [];
    %params.fitc_stringency = 2;
    %params.txred_stringency = 2;
    for i = 1:2:length(param_import)
        switch(param_import{i})
            case {'fitc_min_stretch'}
                params.fitc_min_stretch = str2num(param_import{i+1});
            case {'fitc_min'}
                params.fitc_min = str2num(param_import{i+1});
            case {'fitc_max_stretch'}
                params.fitc_max_stretch = str2num(param_import{i+1});
            case {'fitc_max_foci_size'}
                params.fitc_max_foci_size = str2num(param_import{i+1});
            case {'fitc_k_val'}
                params.fitc_k_val = str2num(param_import{i+1});
            case {'txred_k_val'}
                params.txred_k_val = str2num(param_import{i+1});
            case {'fitc_stringency'}
                params.fitc_stringency = str2num(param_import{i+1});
            case {'txred_stringency'}
                params.txred_stringency = str2num(param_import{i+1});
            case {'txred_min_stretch'}
                params.txred_min_stretch = str2num(param_import{i+1});
            case {'txred_min'}
                params.txred_min = str2num(param_import{i+1});
            case {'txred_max_stretch'}
                params.txred_max_stretch = str2num(param_import{i+1});
            case {'txred_max_foci_size'}
                params.txred_max_foci_size = str2num(param_import{i+1});
            case {'nuc_rad'}
                params.nuc_rad = str2num(param_import{i+1});
            case {'dapi'}
                params.dapi_index = str2num(param_import{i+1});
            case {'fitc'}
                params.fitc_index = str2num(param_import{i+1});
            case {'txred'}
                params.txred_index = str2num(param_import{i+1});
            case {'cy5'}
                params.cy5_index = str2num(param_import{i+1});
            case {'dapi_mask'} %If >0, then by default will use channel dapi_mask as a binary mask for cutting out DAPI field. If string is passed, will assume this string is in dapi_mask name
                params.mask = str2num(param_import{i+1});
            case {'media'}
                params.media = param_import{i+1};
            case {'mode'}
                params.mode = param_import{i+1};
            case {'main_spot'}
                params.spot_index = param_import{i+1};
            case {'measurement'}
                params.msIDs = param_import{i+1};
                temp =textscan(params.msIDs,'%s','delimiter',',');
                params.msIDs = temp{1};
            case {'dxy'}
                params.dxy = str2num(param_import{i+1});
            case {'dz'}
                params.dz = str2num(param_import{i+1});
            otherwise
                  params = feval(@setfield,params,param_import{i},str2num(param_import{i+1}));
        end
    end
catch
    prompt={'Media (i.e. delta,meta,lsm,zvi):','dapi_mask (i.e. Yes, no)','dapi chan: (i.e. 1, 2?}',...
        'fitc chan: (i.e. 1, 2?}','txred chan: (i.e. 1, 2?}','min_stretch','max_stretch','stringency (i.e. 1, 2)',...
        'nuclear radius','min_foci'};
    name='Input for Spot Detection';
    defaultanswer={'meta','no','1','2','','200','4095','2','40','100'};
    answer=inputdlg(prompt,name,1,defaultanswer);
    params.media = answer{1};
    params.mask = strcmp(lower(answer{2}),'yes');
    params.dapi_index = str2num(answer{3});
    params.fitc_index = str2num(answer{4});
    params.txred_index = str2num(answer{5});
    params.fitc_min_stretch = str2num(answer{6});
    params.txred_min_stretch = str2num(answer{6});
    params.fitc_max_stretch = str2num(answer{7});
    params.txred_max_stretch = str2num(answer{7});
    params.fitc_stringency = str2num(answer{8});
    params.txred_stringency = str2num(answer{8});
    params.nuc_rad = str2num(answer{9});
    params.fitc_min = str2num(answer{10});
    params.txred_min = str2num(answer{10});
    params.msIDs = 'size,mean,perc=0.5,rdna';
end

% Independently of input method, fill stringency criteria

if isfield(params, 'fitc_stringency')
    switch (params.fitc_stringency)
        case 1
            params.fitc_k_val = 10;
            params.fitc_res_val = 3;
            params.fitc_max_foci_size = 10;
            params.fitc = 1;    %need this to use previously written output functions
        case 2
            params.fitc_k_val = 2;
            params.fitc_res_val = 3;
            params.fitc_max_foci_size = 10;
            params.fitc = 1;    %need this to use previously written output functions
        case 3
            params.fitc_k_val = 0.5;
            params.fitc_res_val = 3;
            params.fitc_max_foci_size = 5;
            params.fitc = 1;    %need this to use previously written output functions
    end
end
if isfield(params, 'txred_stringency')
    switch (params.txred_stringency)
        case 1
            params.txred_k_val = 10;
            params.txred_res_val = 3;
            params.txred_max_foci_size = 10;
            params.txred = 1;    %need this to use previously written output functions
        case 2
            params.txred_k_val = 2;
            params.txred_res_val = 3;
            params.txred_max_foci_size = 10;
            params.txred = 1;    %need this to use previously written output functions
        case 3
            params.txred_k_val = 0.5;
            params.txred_res_val = 3;
            params.txred_max_foci_size = 5;
            params.txred = 1;    %need this to use previously written output functions
    end
end

%default params for output. may not be needed.
%     params.nuc_out.Filename= 1;
%     params.nuc_out.NucID= 1;
%     params.nuc_out.NucXPosition= 1;
%     params.nuc_out.NucYPosition= 1;
%     params.nuc_out.StageXPosition= 1;
%     params.nuc_out.StageYPosition= 1;
%     params.nuc_out.DapiIntensity= 1;
%     params.nuc_out.DapiStdDev= 1;
%     params.nuc_out.Dose= 1;
%     params.nuc_out.Time= 1;
%     params.nuc_out.NucArea= 1;
%     params.nuc_out.NucVolume= 1;
%     params.nuc_out.FociCount= 1;
%     params.nuc_out.P2A= 1;
%     params.nuc_out.Rdna= 1;
%     params.nuc_out.Rgrad= 1;
%     params.nuc_out.rRdna= 1;
%     params.nuc_out.rRgrad= 1;
%     params.nuc_out.RdnaRatio= 1;
%     params.nuc_out.RgradRatio= 1;
%     params.nuc_out.FitcIntensity= 1;
%     params.nuc_out.FitcStdDev= 1;
%     params.nuc_out.TxRedIntensity= 1;
%     params.nuc_out.TxRedStdDev= 1;
%
%     params.foci_out.Filename= 1;
%     params.foci_out.NucID= 1;
%     params.foci_out.FociID= 1;
%     params.foci_out.Dose= 1;
%     params.foci_out.Time= 1;
%     params.foci_out.FociSize= 1;
%     params.foci_out.FociMeanIntensity= 1;
%     params.foci_out.FociMaxIntensity= 1;
%
%     params.img_out.NucMask= 1;
%     params.img_out.EdgeMask= 1;
%     params.img_out.SpotMask= 1;
%     params.img_out.FociMask= 1;
%     params.img_out.EdgeOverlay= 1;
%     params.img_out.FinalOverlay= 1;

% save 'params.mat' params;
end