function [k, X]=SuperBlockGridCreation(k, xmax, ymax, X, plotit, n_max)
%% SuperBlockGridCreation return the super block grid (k.sb)
%
% INPUT:
%
% * k       : kriging information
% * nx,ny   : number of cell (in x and y)
% * x-ymax  : regular grid max length
% * X       : Primary variable
%
% OUTPUT:
%
% * k       : kriging information with sb grid info
% * X       : Primary variable with assign sb position for all pt
%
% * *Author:* Raphael Nussbaumer (raphael.nussbaumer@unil.ch)
% * *Date:* 02.02.2015


k.sb.dx=xmax/k.sb.nx;
k.sb.x=k.sb.dx/2:k.sb.dx:xmax; % x will start not at 0 but at dx/2
k.sb.dy=ymax/k.sb.ny;
k.sb.y=k.sb.dy/2:k.sb.dy:ymax;


% Creation of the superblock grid windows search
el_max = ceil(max(min([k.covar(1).range*k.wradius; xmax ymax])./[k.sb.dx k.sb.dy]));
[el_X, el_Y] = meshgrid(-el_max:el_max);% grid of searching windows in supergrid unit. this is a quadrant
[el_X_T, el_Y_T]=rotredtrans(el_X*k.sb.dx, el_Y*k.sb.dy, k.covar(1).azimuth, k.covar(1).range); % transforms the grid in unit
el_dist = sqrt(el_X_T.^2 + el_Y_T.^2); % distence from the point 0,0
k.el_X=el_X(el_dist(:)<=k.wradius); 
k.el_Y=el_Y(el_dist(:)<=k.wradius); % All point inside the windows search

% Assignment of all primary data to their super grid
X.sb_x = min([round((X.x-k.sb.x(1))/k.sb.dx +1)'; k.sb.nx*ones(1,X.n)]); % min is to avoid pb when X.x=k.sb.n
X.sb_y = min([round((X.y-k.sb.y(1))/k.sb.dy +1)'; k.sb.ny*ones(1,X.n)]);

% Creation of a mask to select all point which belong to the windows search
% for all superblock grid cell.

k.sb.mask=false(k.sb.ny,k.sb.nx,X.n);

for i=1:k.sb.nx
    for j=1:k.sb.ny % for each cell of the supergrid...
        kt=false(X.n,1);
        for u=1:length(k.el_X) %.. look at all point...
            kt(i+k.el_X(u)==X.sb_x' & j+k.el_Y(u)==X.sb_y')=true;
        end
        k.sb.mask(j,i,kt)=true;
%         
%         [~,id] = sort(sqrt( (X.x(kt)-k.sb.x(i)).^2 + (X.y(kt)-k.sb.y(j)).^2 ));
%         kt2=find(kt);
%         k.sb.mask(j,i,kt2(id(1:min(n_max,numel(id)))))=true;
    end
end


if plotit
    i=5; j=6;
    windows=false(k.sb.ny,k.sb.nx);
    for u=1:length(k.el_X_s) %.. look at all point...
        for q=1:4 %... and assign it to the corresponding quadrant
            if i+k.qs(q,1)*k.el_X_s(u)<=k.sb.nx && j+k.qs(q,2)*k.el_Y_s(u)<=k.sb.ny && i+k.qs(q,1)*k.el_X_s(u)>=1 && j+k.qs(q,2)*k.el_Y_s(u)>=1% check to be inside the grid
                windows(i+k.qs(q,1)*k.el_X_s(u), j+k.qs(q,2)*k.el_Y_s(u))=true;
            end
        end
    end
    
    if plotit
        imagesc(k.sb.x,k.sb.y,windows)
        mesh([0 k.sb.x+k.sb.dx/2],[0 k.sb.y+k.sb.dy/2],zeros(k.sb.nx+1,k.sb.ny+1),'EdgeColor','k','facecolor','none')
        plot(X.x, X.y,'d')
        plot(k.sb.x(i), k.sb.y(j),'or')
        plot(X.x(k.sb.mask(j,i,:)), X.y(k.sb.mask(j,i,:)),'x','lineWidth',3)
        plot([X.x'; k.sb.x(X.sb_x)],[X.y'; k.sb.y(X.sb_y)])
        axis equal tight
        legend('Super Grid','Hard data',['Center of grid ' num2str(i) ';' num2str(j)],['Hard data selected with grid ' num2str(i) ';' num2str(j)])
    end
end

