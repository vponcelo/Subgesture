function [X,Y] = getFullDev(Xdevp,Ydevp)
% return full development (training + validation) sequences
% output:
%   X: Concatenated training + validation data
%   Y: Concatenated training + validation labels
% input
%   Xdevp: Training and validation data split
%   Ydevp: Training and validation labels split

X = toMat(Xdevp);

Y.Lfr = [Ydevp{1}.Lfr Ydevp{2}.Lfr];
Y.L = [Ydevp{1}.L Ydevp{2}.L];
Y.seg = Ydevp{1}.seg;
Y.seg(end) = Ydevp{1}.seg(end) + 1;
segs = Ydevp{2}.seg(2:end) + Y.seg(end)-1;
Y.seg = [Y.seg segs];
