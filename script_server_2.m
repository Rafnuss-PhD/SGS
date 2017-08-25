parm.k.covar.model = '';
parm.k.covar.azimuth = 0;
parm.k.covar.c0 = 1;
parm.k.covar.alpha = 1;
parm.seed_path = Inf;
parm.seed_search = 'shuffle';


parm.k.covar.range0 = [15 15];
parm.k.wradius = Inf;
parm.mg = 0;
nx = 2^6+1;
ny = 2^6+1;
parm.k.nb = 20;


vario = {'exponential', 'gaussian', 'spherical', 'hyperbolic','k-bessel', 'cardinal sine'};

% True
[Y, X] = ndgrid(1:nx,1:ny);
XY = [Y(:) X(:)];

parpool(6)

parfor v=1:numel(vario)
    covar = parm.k.covar;
    covar(1).model = vario{v};
    covar = kriginginitiaite(covar);
    DIST = squareform(pdist(XY*covar.cx));
    CY_true{v} = kron(covar.g(DIST), covar.c0);
end
err_frob_fx = @(CY,v) sqrt(sum((CY(:)-CY_true{v}(:)).^2)) / sum((CY_true{v}(:).^2));


% Simulation
N = 10;
vario_g=cell(numel(vario),1);
D = pdist2(XY,XY);
h = unique(D);
h = h(h<=parm.k.covar.range0(1)*3);

parfor v=1:numel(vario)
    parm1=parm;
    parm1.k.covar.model = vario{v};
    CY=repmat({nan(ny*nx,nx*ny,2)},N,1);
    eta{v}=cell(N,1);
    nn=cell(N,1);
    for n=1:(2^(N-1))
        vec = de2bi(n-1,N);
        CY{1}(:,:,vec(1)+1) = full(SGS_varcovar(nx,ny,parm1));
        eta{v}{1} = [eta{v}{1}; err_frob_fx(CY{1}(:,:,vec(1)+1),v)];
        nn{1} = [nn{1}; 2^(1-1)];
        i=1;
        while (vec(i)==1)
            CY{i+1}(:,:,vec(i+1)+1) = mean(CY{i},3);
            eta{v}{i+1} = [eta{v}{i+1}; err_frob_fx(CY{i+1}(:,:,vec(i+1)+1),v)];
            nn{i+1} = [nn{i+1}; 2^(i)];
            i=i+1;
        end
        disp(['N: ' num2str(n) ])
    end
    
    vario_g{v} = nan(N,numel(h));
    for n=1:N
        for i=1:numel(h)
            id = D ==h(i);
            id2 = CY{n}(id);
            vario_g{v}(n,i) = sum(id2)/sum(id(:));
        end
    end
    %save(['frobenium_',vario{v},'_20_512_MG'])
end


nn=cell(N,1);
for n=1:(2^(N-1))
    vec = de2bi(n-1,N);
    nn{1} = [nn{1}; 2^(1-1)];
    i=1;
    while (vec(i)==1)
        nn{i+1} = [nn{i+1}; 2^(i)];
        i=i+1;
    end
end

save(['./cst_path_paper/frobenium_65n_20k_512N_0MG'],'eta','vario','vario_g','nn','parm','nx','ny','h');
