clear all;

parm.k.covar.model = 'spherical';
parm.k.covar.range0 = [10 10];
parm.k.covar.azimuth = 0;
parm.k.covar.c0 = 1;
parm.k.covar.alpha = 1;
parm.saveit = false;
parm.n_real  = 1;
parm.seed_path = 'shuffle';
parm.seed_search = 'shuffle';
parm.seed_U = 'shuffle';
parm.varcovar = 0;
parm.k.wradius = 2;
parm.k.lookup = true;

parm.mg = 0;

nx = 50; % no multigrid
ny = 50;
parm.k.nb = 4;


% Run the program
% parpool(4);
% profile on
tic
Res = SGS_cst_par(nx,ny,parm);
Res = (Res-mean(Res(:)))./std(Res(:));
toc
%profile viewer

% Plot
figure(1); clf;
imagesc(Res);
xlabel('x'); ylabel('y'); colorbar; axis equal



%% Conditional

parm.mg = 0;
parm.k.nb = 4;
parm.k.lookup = true;

% Conditional
grid.x=1:nx;
grid.y=1:ny;
addpath('./data_gen/')
field=fftma_perso(parm.k.covar, grid);
prim = sampling_pt(grid,field,2,0);


Res = SGS_cst_par_cond(nx,ny,prim,parm);
Res = (Res-mean(Res(:)))./std(Res(:));

% Plot
figure(1); clf;
subplot(1,2,1); hold on; imagesc(field);scatter(prim.x,prim.y,[],prim.d,'filled','MarkerEdgeColor','k'); colorbar; axis tight equal; xlabel('x'); ylabel('y'); 
subplot(1,2,2); hold on; imagesc(Res); scatter(prim.x,prim.y,[],prim.d,'filled','MarkerEdgeColor','k'); colorbar; axis tight equal; xlabel('x'); ylabel('y'); 

figure(2); clf; hold on;
plot(variogram_gridded_perso(Res))
plot(variogram_gridded_perso(field))
covar=kriginginitiaite(parm.k.covar);
plot(1-covar.g((0:100)/covar.range(1)),'--k')
xlim([1 covar.range(1)*3])
