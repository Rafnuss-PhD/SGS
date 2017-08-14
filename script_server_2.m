addpath('./functions/.')


% General setting
Z.d=[];Z.x=[];Z.y=[];Z.n=0; % no hard data

parm.k.covar(1).model = 'spherical';
parm.k.covar(1).range0 = [15 15];
parm.k.covar(1).azimuth = 0;
parm.k.covar(1).c0 = 1;
parm.k.covar(1).alpha = 1;
% parm.k.covar(2).model = 'nugget';
% parm.k.covar(2).range0 = [0 0];
% parm.k.covar(2).azimuth = 0;
% parm.k.covar(2).c0 = 0.02;

parm.saveit = 0;
parm.nscore = 0;
parm.seed = 'shuffle';
parm.k.method = 'minkmex';
parm.k.nb = [0 0 0 0 0; 20 0 0 0 0];
parm.g.scale=[6;6];

parm.path = 'linear';
parm.path_random = 1;
parm.varcovar = 1;
parm.par = 0;


vario = {'exponential', 'gaussian', 'spherical', 'hyperbolic','k-bessel', 'cardinal sine'};


% True
grid_gen.x = 0:2^parm.g.scale(1);
grid_gen.y = 0:2^parm.g.scale(2);
[grid_gen.X, grid_gen.Y] = meshgrid(grid_gen.x,grid_gen.y);
XY = [grid_gen.X(:) grid_gen.Y(:)];
D = pdist2(XY,XY);
h = unique(D);
h=h(h<parm.k.covar(1).range0(1)*2);

for v=1:numel(vario)
    parm1=parm;
    parm1.k.covar(1).model = vario{v};
    covar = kriginginitiaite(parm1.k.covar);
    vario_true{v}=covar(1).g(h./covar.range(1));
    %plot(h,covar(1).g(h./covar.range(1)))
end




% Simulation
N=512;
parpool(48);
for v=1:numel(vario)
    parm1=parm;
    parm1.k.covar(1).model = vario{v};
    vario_sim1 = nan(N,numel(h));
    parfor i=1:N
        Res = SGS(Z,parm1);
        CY = (Res{end}.lambda) \ transpose(inv(Res{end}.lambda));
        vario_sim2 = nan(1,numel(h));
        for j=1:numel(h)
            id = D==h(j);
            id2 = CY(id);
            vario_sim2(j) = sum(id2)/sum(id(:));
        end
        vario_sim1(i,:) = vario_sim2;
    end
    vario_sim{v} = vario_sim1;
    
    disp(['Vario done: ' num2str(v)])
end

save('vario_all','vario_true','vario_sim','parm','vario','N','h')
