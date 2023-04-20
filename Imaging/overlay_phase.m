function ov_img = overlay_phase(phase,gray,color)
% ov_img = overlay_phase(phase,gray,color)
% Will overlay gray image as a color image over the phase image.
% color can be blue, red or green.
% Default: Green
%
% February 2007, Sylvain Costes, LBNL

if ~exist('color')
    color = 'green';
end

phase = dip_image(uint8(stretch(phase,0,100,0,127)));
gray = dip_image(uint8(stretch(gray,0,100,0,127)));

switch lower(color)
    case 'red'
        ov_img = newimar(phase+gray,phase,phase);
    case 'blue'
        ov_img = newimar(phase,phase,phase+gray);
    otherwise
        ov_img = newimar(phase,phase+gray,phase);
end

ov_img = colorspace(ov_img,'rgb');