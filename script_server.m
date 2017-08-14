addpath('./functions/.')

% General setting
Z.d=[];Z.x=[];Z.y=[];Z.n=0; % no hard data

parm.k.covar(1).model = 'spherical';
parm.k.covar(1).range0 = [15 15];
parm.k.covar(1).azimuth = 0;
parm.k.covar(1).c0 = .98;
parm.k.covar(1).alpha = 1;
parm.k.covar(2).model = 'nugget';
parm.k.covar(2).range0 = [0 0];
parm.k.covar(2).azimuth = 0;
parm.k.covar(2).c0 = 0.02;

parm.saveit = 0;
parm.nscore = 0;
parm.par = 0;
parm.seed = 'default';



%% Computational Saving cst path
parm.g.scale=[7;7];

parm.par = 0;

parm.k.wradius = 2;
parm.k.method = 'sbss';
parm.k.nb = [0 0 0 0 0; 30 0 0 0 0];

parm.n_realisation  = 1;
parm.cstk = 1;

[Res,t] = SGS(Z,parm);

t
sum(t.krig)/t.global
save('./cst_path_paper/T_12x12n_30k_1m_MINK','parm','t')
