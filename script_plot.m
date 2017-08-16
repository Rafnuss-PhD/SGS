clear all;

parm.k.covar(1).model = 'spherical';
parm.k.covar(1).range0 = [15 15];
parm.k.covar(1).azimuth = 0;
parm.k.covar(1).c0 = 1;
parm.k.covar(1).alpha = 1;
% parm.k.covar(2).model = 'nugget';
% parm.k.covar(2).range0 = [0 0];
% parm.k.covar(2).azimuth = 0;
% parm.k.covar(2).c0 = 0.02;
parm.saveit = false;
parm.n_real  = 1;
parm.seed_path = 'default';
parm.seed_search = 'shuffle';
parm.seed_U = 'default';
parm.varcovar = 0;
parm.k.wradius = 10;

parm.mg = 1;

nx = 4; % no multigrid
ny = 4;
parm.k.nb = 1;

% Run the program
% parpool(4);
% profile on
tic
Res = SGS_par(nx,ny,parm);
toc
%profile viewer

% Plot
figure(1); clf;
imagesc(Res);
xlabel('x'); ylabel('y'); colorbar;



