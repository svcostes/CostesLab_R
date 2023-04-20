function str = int2str_format(i,ndigit)
% int2str will turn an integer into a formated string.
% input: i,ndigit
% i: integer to change
% dec: # of digit to the left of decimal point wanted. 
% Example
% >> int2str_format(12,4)
% ans =
% 0012
% Sylvain Costes, LBNL, April 2005

r = i;
str = '';
for k = 1:ndigit
    digit(k) = fix(r/10^(ndigit-k));
    str = sprintf('%s%d',str,digit(k));
    r = r - digit(k)*10^(ndigit-k);
end

