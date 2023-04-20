function [img_array,metadata,timeLapse] = bfopen3(id)

% A script for opening microscopy images in MATLAB using Bio-Formats.
%
% Portions of this code were adapted from:
% http://www.mathworks.com/support/solutions/data/1-2WPAYR.html?solution=1-2WPAYR
%
% Excluding figure creation, this method is ~1.5x-2.5x slower than
% Bio-Formats's command line showinf tool (MATLAB R14 vs. java 1.6.0_03),
% due to overhead from reshaping arrays and converting pixel types.
%
% Thanks to Ville Rantanen for his performance improvements and ideas.
% Thanks to Brett Shoelson of The MathWorks for his excellent suggestions.
%
% To install, down load loci_tools.jar from:
%   http://www.loci.wisc.edu/ome/formats.html
% Internet Explorer sometimes erroneously renames the Bio-Formats library
% to loci_tools.zip. If this happens, rename it back to loci_tools.jar.
% Place loci_tools.jar and this script (bfopen.m) in your MATLAB work folder.
%
% Modified August 2009 for DIPimage. Output will return an dipimage array.
% Each image in the array is for another channel. The image is
% multidimension: (W,H,D,T): Width: W, Height: H, Z depth: D, time length:
% T.
%
% Modified April 2011 for supporting the metadata time lapse information
% returns an array of time lapse for a time series image. (By Charlie)
%
% Sylvain Costes, Lawrence Berkeley Lab.
%


% load Bio-Formats library into MATLAB environment
if isunix  % Mac & Linux
    [idum,hostname] = system('hostname');
    hostname=hostname(1:end-1); % truncate the newline
    OS=computer; % figure out OS type
    if strcmp(OS, 'MACI64')  % if it's Mac
        if or(strcmp(hostname, 'orona'),strcmp(hostname, 'dumbledore'))
            javaaddpath('/Users/genuser/PrivateServers/NSCOR/MyMatlab/loci_tools.jar'); % For Mac on orona
        else
            javaaddpath('/Volumes/MyMatlab/loci_tools.jar'); % For Mac
        end
    else  % Linux
        javaaddpath('/NSCOR_MyMatlab/loci_tools.jar') % For linux
    end
else  % Windows
    javaaddpath('\\radbio\data2\NSCOR_MyMatlab\loci_tools.jar');
end
% Alternately, you can add the library to MATLAB's static class path:
%   1. Type "edit classpath.txt" at the MATLAB prompt.
%   2. Go to the end of the file, and add the path to your JAR file
%      (e.g., C:/Program Files/MATLAB/work/loci_tools.jar).
%   3. Save the file and restart MATLAB.

% initialize logging
loci.common.DebugTools.enableLogging('INFO');

r = loci.formats.ChannelFiller();
r = loci.formats.ChannelSeparator(r);
% comment the following line if you are using a stable release
r = loci.formats.gui.BufferedImageReader(r);

% uncomment the following line to enable grouping of similarly
% named files into a single dataset based on file numbering
% r = loci.formats.FileStitcher(r);

tic
r.setId(id);
numSeries = r.getSeriesCount();
result = cell(numSeries, 2);
metadata = r.getMetadata;

% temporary solution for weirdness of Thambi's .lsm file
% added by Charlie
%[head, ext] = strread(id, '%s%s', 'delimiter', '.');
% if strcmp(ext, 'lsm')
%    numSeries = 2; 
% end

for s = 1:numSeries
    fprintf('Reading series #%d', s);
    r.setSeries(s - 1);
    w = r.getSizeX();
    h = r.getSizeY();
    shape = [w h];
    numImages = r.getImageCount();
    imageList = cell(numImages, 2);
    sizeZ = r.getSizeZ();
    sizeC = r.getSizeC();
    sizeT = r.getSizeT();
    
    img_array = newimar(sizeC);
    for i = 1:sizeC
        img_array{i} = newim(w,h,sizeZ,sizeT);
    end
    
    for i = 1:numImages
        fprintf('.');
        img = r.openImage(i - 1);
        % convert Java BufferedImage to MATLAB image
        pix = img.getData.getPixels(0, 0, w, h, []);
        arr = reshape(pix, shape)';
        % build an informative title for our figure
        label = id;
        if numSeries > 1
            qs = int2str(s);
            label = [label, '; series ', qs, '/', int2str(numSeries)];
        end
        if numImages > 1
            qi = int2str(i);
            label = [label, '; plane ', qi, '/', int2str(numImages)];
            if r.isOrderCertain()
                lz = 'Z';
                lc = 'C';
                lt = 'T';
            else
                lz = 'Z?';
                lc = 'C?';
                lt = 'T?';
            end
            zct = r.getZCTCoords(i - 1);
            if sizeZ > 1
                qz = int2str(zct(1) + 1);
                label = [label, '; ', lz, '=', qz, '/', int2str(sizeZ)];
            end
            if sizeC > 1
                qc = int2str(zct(2) + 1);
                label = [label, '; ', lc, '=', qc, '/', int2str(sizeC)];
            end
            if sizeT > 1
                qt = int2str(zct(3) + 1);
                label = [label, '; ', lt, '=', qt, '/', int2str(sizeT)];
            end
            img_array{zct(2)+1}(:,:,zct(1),zct(3)) = dip_image(arr);
        else   % fixed dip_image not created when numImages == 1
            img_array{i} = dip_image(arr);
        end
        % save image plane and label into the list
        imageList{i, 1} = arr;
        imageList{i, 2} = label;
    end
    for i = 1:sizeC
        img_array{i} = squeeze(img_array{i});
    end
    % extract metadata table for this series
    metadataList = r.getMetadata();
    % get the time lapse information, if possible
    if sizeT > 1
        CONVERSION_CONST = 1.5912e+06; % Based on reverse engineering from AxioVision
        incrementFactor = sizeC*sizeZ;
        endNum = sizeC*sizeZ*sizeT;
        timeLapse = [];
        for time=1:incrementFactor:endNum
            tmp= str2double(metadataList.get(['Timestamp ' num2str(time)]));
            timeLapse = [timeLapse tmp];
        end
    else
        timeLapse = [];
    end
    % save images and metadata into our master series list
    result{s, 1} = imageList;
    result{s, 2} = metadataList;
    result{s, 3} = timeLapse;
    fprintf('\n');
end
% force java garbage collection after the reading, to save some heap space
java.lang.System.gc
r.close;
toc

