function pt = sampling_pt(grid,field,method,k)
% SAMPLING_PT simulate the measurement of high resolution point of the a
% field. Two method are implemented (borehole or random)
% INPUT:
%       - grid:     grid of the matrix to generate (see Import_dat for more info)
%       - field:   true g field
%       - method:   choice of method : 1. borehole | 2. random generated
%       - plotit:   1 or 0 to disply a plot or not
%
% OUTPUT:
%       - pt:        points measurement (pt.x, pt.y, pt.d)
%
% Author: Raphael Nussbaumer
% date : January 2014
% need to do : add assert() for input, flexible number of input var

if isempty(field)
    warning('field is empty')
    pt=[];
    return
end

switch method
    case 1 % select all point in k boreholes equaly spaced
        
        %  positions of conditioning data
        [y_id, x_id] = ndgrid(1:numel(grid.y), round(linspace(1,numel(grid.x),k)));
        pt.x = grid.x(x_id(:))';
        pt.y = grid.y(y_id(:))';
        
        % Create the input data
        pt.id = sub2ind(size(field),y_id(:),x_id(:));
        pt.d = field(pt.id);
        

    case 2 % select k random point on the mesh
        % rng(123456)
        [pt.d,pt.id] = datasample(field(:),k,'Replace',false);
        [pt_j,pt_i] = ind2sub(size(field),pt.id);
        pt.x=grid.x(pt_i);
        pt.y=grid.y(pt_j);

end
pt.n=numel(pt.d);
end
