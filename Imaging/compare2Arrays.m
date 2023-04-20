% Function to compare two arrays and return the common elements indices
% INPUT:
% A & B are two vectors
% OUTPUT:
% X & Y are two vectors such that A(X) == B(Y)
%EXAMPLE:
% A = [ 20 10 100 120 80 50 90 71 85 1201 ];
% B=randperm(100);
% [X Y] = compare2Arrays(A,B);
% max(A(X)-B(Y))
%
% Copyrights:
% Author: Chenna Krishna Varri
% : Computer Vision & Robotics Laboratory, George Mason University.
%Modified: July 2006

function [X Y] = compare2Arrays(A,B)

% Transposing the arrays if they are not comlumn vectors
[m,n]=size(A);
if n>=m
A=A';
end
[m,n]=size(B);
if n>=m
B=B';
end

C=repmat(A,1,size(B,1));
D=C-repmat(B',size(A,1),1);
[X,Y]=find(D==0);
