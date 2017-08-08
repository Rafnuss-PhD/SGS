function [x,y]=rotredtrans(x,y,ang,range)
% function [cx,rot]=trans(cx,model,im);
%
% TRANS is called from COKRI2. It takes as input original coordinates and
%       return the rotated and reduced coordinates following specifications
%       described in model(im,:)
% Rotations are all performed anticlockwise with the observer located on the positive side of
% the axis and looking toward the origin. In 3D, rotations are performed first along z,
% then along rotated y and then along twice rotated x.
% Author: D. Marcotte
% Version 2.1  97/aug/18
[i,j]=size(x);
assert(all([i,j]==size(y)),'need to be same size')

% perform rotation counterclockwise
cang=cos(ang/180*pi); sang=sin(ang/180*pi);
rot=[cang,-sang;sang,cang];

% rotation is performed around z, y and x in that order, the other coordinates are left unchanged.
% perform contractions or dilatations (reduced h)

t=[x(:) y(:)]*rot/diag(range);

x=reshape(t(:,1),i,j);
y=reshape(t(:,2),i,j);
