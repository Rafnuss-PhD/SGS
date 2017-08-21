%% SECTION: COMPUTATIONAL SAVING

clear all;

% General setting
parm.k.covar(1).model = 'spherical';
parm.k.covar(1).azimuth = 0;
parm.k.covar(1).c0 = 1;
parm.k.covar(1).alpha = 1;
parm.seed_path = 'default';
parm.seed_search = 'shuffle';
parm.seed_U = 'default';

parm.k.covar(1).range0 = [15 15];
parm.k.wradius = 1;
parm.n_real=1;
parm.saveit = 0;

parm.mg = 1;


N=[255 511 1023 2047 4095 8191];
K=[20 52 108];

nx = N(1); % no multigrid
ny = N(1);
parm.k.nb = K(1);

[~,t] = SGS_cst_par(nx,ny,parm)

% [~,t] = SGS_trad(nx,ny,parm)


% Run the program



save(['./cst_path_paper/T_' num2str(nx) 'N_' num2str(parm.k.nb) 'K' ],'parm','t')
