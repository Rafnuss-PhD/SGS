%%
% General setting
parm.k.covar.model = 'spherical';
parm.k.covar.azimuth = 0;
parm.k.covar.c0 = 1;
parm.k.covar.alpha = 1;
parm.seed_path = 9;
disp(num2str(parm.seed_path))
parm.seed_search = 'default';


parm.k.covar.range0 = [15 15];
parm.k.wradius = 1;
parm.mg = 1;
nx = 2^7+1; % no multigrid
ny = 2^7+1;
parm.k.nb = 20;


% True
[Y, X]=ndgrid(1:nx,1:ny);
XY = [Y(:) X(:)];
covar = kriginginitiaite(parm.k.covar);
DIST = squareform(pdist(XY*covar.cx));
CY_true = kron(covar.g(DIST), covar.c0);
err_frob_fx = @(CY) sqrt(sum((CY(:)-CY_true(:)).^2)) / sum((CY_true(:).^2));


N=5;
CY=repmat({nan(ny*nx,nx*ny,2)},N,1);
eta=cell(N,1);
nn=cell(N,1);
for n=1:(2^(N-1))
    vec = de2bi(n-1,N);
    CY{1}(:,:,vec(1)+1) = full(SGS_varcovar(nx,ny,parm));
    eta{1} = [eta{1}; err_frob_fx(CY{1}(:,:,vec(1)+1))];
    nn{1} = [nn{1}; 2^(1-1)];
    i=1;
    while (vec(i)==1)
        CY{i+1}(:,:,vec(i+1)+1) = mean(CY{i},3);
        eta{i+1} = [eta{i+1}; err_frob_fx(CY{i+1}(:,:,vec(i+1)+1))];
        nn{i+1} = [nn{i+1}; 2^(i)];
        i=i+1;
    end
    disp(['N: ' num2str(n) ])
end

save(['./cst_path_paper/frobenium_129n_20k_16N_MG_' num2str(parm.seed_path) 'CSTP'],'eta','nn','parm','nx','ny');




