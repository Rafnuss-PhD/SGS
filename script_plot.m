
% Add the function path
addpath(genpath('./functions/.'))


% Define the general setting
Z.d=[];Z.x=[];Z.y=[];Z.n=0; % no hard data
parm.g.scale=[6 ;6]; % no multigrid
parm.g.max=[100 100];
parm.k.covar(1).model = 'spherical';
parm.k.covar(1).range0 = [15 15];
parm.k.covar(1).azimuth = 0;
parm.k.covar(1).c0 = .98;
parm.k.covar(1).alpha = 1;
parm.k.covar(2).model = 'nugget';
parm.k.covar(2).range0 = [0 0];
parm.k.covar(2).azimuth = 0;
parm.k.covar(2).c0 = 0.02;
parm.saveit = false;
parm.nscore = 0;
parm.par = 0;
parm.n_realisation  = 1;
parm.cstk = true;
parm.seed = 'default';
parm.varcovar = 0;

parm.k.wradius=Inf;

% Run the program
[Res,t] = SGS(Z,parm);


% Plot
figure(1); clf;
imagesc(Res{1}.x, Res{1}.y, Res{end}.m{1});
xlabel('x'); ylabel('y'); colorbar;
