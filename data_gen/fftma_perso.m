function field_f=fftma_perso(covar, grid)

%% Create super grid 
grid_s.x_min = grid.x(1);
grid_s.x_max = grid.x(end)*3;
grid_s.y_min = grid.y(1);
grid_s.y_max = grid.y(end)*3;

if ~isfield(grid, 'dx')
    grid_s.dx    = grid.x(2)-grid.x(1);
    grid_s.dy    = grid.y(2)-grid.y(1);
    
else
    grid_s.dx    = grid.dx;
    grid_s.dy    = grid.dy;
end

if ~isfield(grid, 'nx')
    grid.nx = numel(grid.x);
    grid.ny = numel(grid.y);
end


covar=kriginginitiaite(covar);


%% Generate field

% addpath('C:\Users\rafnu\Documents\MATLAB\mGstat')
% addpath('C:\Users\rafnu\Documents\MATLAB\mGstat\misc')
% Va=[num2str(gen.covar.c(2)),' Nug(0) + ', num2str(gen.covar.c(1)),' Sph(', num2str(gen.covar.modele(1,2)), ',90,', num2str(gen.covar.modele(1,3)/gen.covar.modele(1,2)) ,')']; % V = �sill Sph(range,rotation,anisotropy_factor)�
% field_g=fft_ma_2d(grid_s.x,grid_s.y,Va);
field_g       = fftma(grid_s.x_min,grid_s.dx,grid_s.x_max,grid_s.y_min,grid_s.dy,grid_s.y_max,covar);


%% Resample the field to initial size
field_p=field_g(grid.ny+1:2*grid.ny,grid.nx+1:2*grid.nx);


%% Adjust the field
field_f = (field_p-mean(field_p(:)))./std(field_p(:));



%% Plot

% figure;imagesc(field_f);axis equal;colorbar;
% figure; hist(field_f(:));
% 
% myfun = @(x,h) semivariogram1D(h,1,x,'sph',0);
% 
% [gamma_x, gamma_y] = variogram_gridded_perso(field_f);
% figure; subplot(1,2,1); hold on;
% id= grid.x<covar.modele(1,2);
% plot(grid.x(id),gamma_x(id),'linewidth',2)
% plot(grid.x(id),myfun(covar.modele(1,2),grid.x(id)),'linewidth',2)
% subplot(1,2,2);hold on
% id= grid.y<covar.modele(1,3);
% plot(grid.y(id),gamma_y(id),'linewidth',2)
% plot(grid.y(id),myfun(covar.modele(1,3),grid.y(id)),'linewidth',2)

end

function [zs]=fftma(xmin,dx,xmax,ymin,dy,ymax,covar)
% [zs]=fftma(xmin,dx,xmax,ymin,dy,ymax,covar)
  
% Copyright (C) 2005 Erwan Gloaguen, Bernard Giroux
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
% 

small=1e-6;

Nx=2*(1+length(xmin:dx:xmax+small));
Ny=2*(1+length(ymin:dy:ymax+small));

Nx2 = Nx/2;
Ny2 = Ny/2;

x = dx*(0:Nx2-1);
x = [x fliplr(-x)]';
y = dy*(0:Ny2-1);
y = [y fliplr(-y)]';

x = kron(ones(Ny,1), x);
y = kron(y, ones(Nx,1));
d = covardm_perso([x y],[0 0],covar);

K = reshape(d,Nx,Ny)';

%%%%%%%%%%%%%%%%%%%%%%%%%%
%On s'assure que la covariance tombe bien à zéro
mk=0;
if min(K(:,1))>1e-6
  Ny=2*Ny;
  mk=1;
end
if min(K(1,:))>1e-6
  Nx=2*Nx;
  mk=1;
end  

if mk==1;
  Nx2 = Nx/2;
  Ny2 = Ny/2;

  x = dx*(0:Nx2-1);
  x = [x fliplr(-x)]';
  y = dy*(0:Ny2-1);
  y = [y fliplr(-y)]';

  x = kron(ones(Ny,1), x);
  y = kron(y, ones(Nx,1));
  d = covardm_perso([x y],[0 0],covar);

  K = reshape(d,Nx,Ny)';
end

% Calcul de G

G=fft2(K).^0.5;

U=fft2(randn(size(K)));

GU=G.*U;

% Transformation de Fourier inverse donnant g*u et z
Z=real(ifft2(GU));
%reconstruction de la simulation sur la taille du champ
if mk==0
  zs=Z((Ny+2)/2+1:end,(Nx+2)/2+1:end);
elseif mk==1
  zs=Z((Ny+2)/2+1:(Ny+2)/2+length(ymin:dy:ymax+small),(Nx+2)/2+1:(Nx+2)/2+length(xmin:dx:xmax+small));
end
end
