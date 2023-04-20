function uiopen(type,direct)
%UIOPEN Present file selection dialog with appropriate file filters.
%
%   UIOPEN presents a file selection dialog.  The user can either choose a
%   file to open or click cancel.  No further action is taken if the user
%   clicks on cancel.  Otherwise the OPEN command is evaluated in the base
%   workspace with the user specified filename.
%
%   These are the file filters presented using UIOPEN.
%
%   1st input argument   Filter List
%   <no input args>      *.m, *.fig, *.mat,
%                        *.mdl         (if Simulink is installed),
%                        *.cdr         (if Stateflow is installed),
%                        *.rtw, *.tmf, *.tlc, *.c, *.h, *.ads, *.adb
%                                      (if Real-Time Workshop is installed),
%                      *.*
%   MATLAB               *.m, *.fig, *.*
%   LOAD                 *.mat, *.*
%   FIGURE               *.fig, *.*
%   SIMULINK             *.mdl, *.*
%   EDITOR               *.m, *.mdl, *.cdr, *.rtw, *.tmf, *.tlc, *.c, *.h, *.ads, *.adb, *.*
%
%   If the first input argument is unrecognized, it is treated as a file
%   filter and passed directly to the UIGETFILE command.
%
%   If the second input argument is true, the first input argument is
%   treated as a filename.
%
%   See also UIGETFILE, UIPUTFILE, OPEN, UILOAD, UIIMPORT.
%
%   Modified by James Chen and Sylvain Costes, LBNL, March 2007 to add stk
%   tif and ics files.
%
%   Copyright 1984-2004 The MathWorks, Inc.
%   $Revision: 1.23.4.7 $  $Date: 2004/10/22 19:45:12 $

if nargin == 0
    type = '';
end

slex=0; sfex=0; rtwex=0;

allML = {'*.m;*.fig;*.mat', 'All MATLAB Files'};
if exist('toolbox/simulink/simulink','dir')
    allML(1)=strcat(allML(1), ';*.mdl');
    slex=1;
end
if exist('toolbox/stateflow/stateflow','dir')
    allML(1)=strcat(allML(1),';*.cdr');
    sfex=1;
end
if exist('toolbox/rtw/rtw','dir')
    allML(1)=strcat(allML(1),';*.rtw;*.tmf;*.tlc;*.c;*.h;*.ads;*.adb');
    rtwex=1;
end

if isempty(type)
    filterList = [
        allML; ...
        {'*.m',   'M-files (*.m)'; ...
        '*.fig', 'Figures (*.fig)'; ...
        '*.mat', 'MAT-files (*.mat)'}
        ];
    if slex
        rows=size(filterList, 1);
        filterList {rows+1, 1}='*.mdl';
        filterList {rows+1, 2}='Models (*.mdl)';
    end
    if sfex
        rows=size(filterList, 1);
        filterList {rows+1, 1}='*.cdr';
        filterList {rows+1, 2}='Stateflow files (*.cdr)';
    end
    if rtwex
        rows=size(filterList, 1);
        filterList {rows+1, 1}='*.rtw;*.tmf;*.tlc;*.c;*.h;*.adb;*.ads';
        filterList {rows+1, 2}='Real-Time Workshop files (*.rtw,*.tmf,*.tlc,*.c,*.h,*.adb,*.ads)';
    end
    rows=size(filterList, 1);
    filterList {rows+1, 1}='*.*';
    filterList {rows+1, 2}='All Files (*.*)';
else
    switch(lower(type))
        case 'matlab'
            filterList = [
                allML; ...
                {'*.m',   'M-files (*.m)'; ...
                '*.fig', 'Figures (*.fig)'}
                ];
            rows=size(filterList, 1);
            filterList {rows+1, 1}='*.*';
            filterList {rows+1, 2}='All Files (*.*)';
        case 'load'
            filterList = [
                {'*.mat', 'MAT-files (*.mat)'}; ...
                allML; ...
                {'*.*',   'All Files (*.*)'}
                ];
        case 'figure'
            filterList = [
                {'*.fig', 'Figures (*.fig)'}; ...
                allML; ...
                {'*.*',   'All Files (*.*)'}
                ];
        case 'simulink'
            filterList = [
                {'*.mdl', 'Model (*.mdl)'}; ...
                allML; ...
                {'*.*',   'All Files (*.*)'}
                ];
        case 'editor'
            
            
            
            MLText = {'*.m', 'All MATLAB Files'};
            MLText(1)=strcat(MLText(1),';*.mdl');
            MLText(1)=strcat(MLText(1),';*.cdr');
            MLText(1)=strcat(MLText(1),';*.rtw;*.tmf;*.tlc;*.c;*.h;*.ads;*.adb');
            
            filterList = [
                MLText;...
                {'*.*',   'All Files (*.*)'}
                ];
            
        otherwise % Add here tif, ics and stk files
            if ~isempty(strfind(lower(type), '.stk'))
                beg = strfind(lower(type),'\');
                en = strfind(lower(type),'.stk');
                temp_img = read_gray_stack(type)
                
            else
                if ~isempty(strfind(lower(type), '.tif'))
                    beg = strfind(lower(type),'\');
                    en = strfind(lower(type),'.tif');
                    temp_img = dip_image(imread(type))
                else
                    if ~isempty(strfind(lower(type), '.ics'))
                        beg = strfind(lower(type),'\');
                        en = strfind(lower(type),'.ics');
                        temp_img = readim(type)
                    else
                        if ~isempty(strfind(lower(type), '.dv'))
                            beg = strfind(lower(type),'\');
                            en = strfind(lower(type),'.dv');
                            temp_img = bfopen3(type)
                        else
                            if ~isempty(strfind(lower(type), '.zvi'))
                                beg = strfind(lower(type),'\');
                                en = strfind(lower(type),'.zvi');
                                temp_img = bfopen3(type)
                            else
                                filterList = type;
                            end
                        end
                    end
                end
            end
            % If image was read, put it into base workspace.
            if exist('temp_img')
                try
                    assignin('base',type(beg(end)+1:en(end)-1),temp_img);
                catch
                    assignin('base','imported_img',temp_img);
                end
                return;
            end
    end
end

if nargin < 2
    direct = false;
end

if direct
    fn = type;
else
    % Is it a .APP or .KEY directory on the Mac?
    % If so, open it properly.
    if strncmp(computer,'MAC',3) && ~iscell(filterList)...
            && isdir(filterList)
        [unused, unused2, ext] = fileparts(filterList);
        if strcmpi(ext, '.app') || strcmpi(ext, '.key')
            unix(['open "' filterList '" &']);
            return;
        end
    end
    [fn,pn] = uigetfile(filterList,'Open');
    if isequal(fn,0)
        return;
    end
    fn = fullfile(pn,fn);
end

try
    % send open requests from editor back to editor
    if strcmpi(type,'editor')
        edit(fn);
    else
        % Is it a MAT-file?
        [unused, unused2, ext] = fileparts(fn);
        if strcmpi(ext, '.mat')
            quotedFile = ['''' strrep(fn, '''', '''''') ''''];
            evalin('caller', ['load(' quotedFile ');']);
            setStatusBar(~isempty(whos('-file', fn)));
            return;
        end
        
        % Look to see if it's an HDF file  If so, don't try to handle it;
        % Pass it off to tools that know what to do.
        out = [];
        fid = fopen(fn);
        if fid ~= -1
            out = fread(fid, 4);
            fclose(fid);
        end
        if length(out) == 4 && sum(out == [14; 3; 19; 1]) == 4
            % The HDF tool is not yet available on GLNXA64.
            if ~strcmp(computer, 'GLNXA64')
                hdftool(fn);
            else
                errordlg('HDF files are not supported on this platform.');
            end
        else
            ans = [];
            % If open creates variables, ans will get populated in here.
            % We need to assign it in the calling workspace later
            open(fn);
            
            if ~isempty(ans)
                vars = ans;
                % Shove the resulting variables into the calling workspace
                status = true;
                if isstruct(vars)
                    names = fieldnames(vars);
                    status = length(names) > 0;
                    for i = 1:length(names)
                        assignin('caller',names{i}, vars.(names{i}));
                    end
                else
                    assignin('caller','ans',vars);
                end
                setStatusBar(status);
            end
        end
    end
catch
    errordlg(lasterr);
end

function setStatusBar(varsCreated)

if varsCreated
    message = 'Variables created in current workspace.';
else
    message = 'No variables created in current workspace.';
end

% The following class reference is undocumented and
% unsupported, and may change at any time.
dt = javaMethod('getInstance', 'com.mathworks.mde.desk.MLDesktop');
if dt.hasMainFrame
    dt.setStatusText(message);
else
    disp(message);
end

