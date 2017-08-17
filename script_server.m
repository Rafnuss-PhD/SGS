clear all;

% General setting
parm.k.covar.model = '';
parm.k.covar.azimuth = 0;
parm.k.covar.c0 = 1;
parm.k.covar.alpha = 1;
parm.seed_path = 'default';
parm.seed_search = 'default';
parm.seed_U = 'default';

parm.k.covar.range0 = [15 15];
parm.k.wradius = 1.3;
parm.mg = 0;
nx = 100; % no multigrid
ny = 100;
parm.k.nb = 30;


vario = {'exponential', 'gaussian', 'spherical', 'hyperbolic','k-bessel', 'cardinal sine'};

% True
[Y, X]=ndgrid(1:nx,1:ny);
XY = [Y(:) X(:)];
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
parfor v=1:numel(vario)
    parm1=parm;
    parm1.k.covar.model = vario{v};
    CY=repmat({nan(ny*nx,ny*nx,2)},N,1);
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
    end
    disp(['Vario: ' num2str(v) '/' num2str(numel(vario))])
end

save(['./cst_path_paper/frobenium_100n_30k_512N'],'eta','vario','nn','parm','nx','ny');

boxplot(cell2mat(eta{v}),cell2mat(nn),'Orientation','horizontal');


%figure
load('result-SGSIM/Constant Path/frobenium_v20_x512_MG');
load('result-SGSIM/Constant Path/frobenium_v20_x512');

figure(2);clf; hold on;
lim=4;
for v=1:numel(vario)
    subplot(numel(vario)/2,2,v); hold on;
    boxplot(cell2mat(eta{v}(1:lim)),cell2mat(nn(1:lim)),'Orientation','horizontal','Color',[0 113 188]/255);
    axis tight;
    xlabel(['Number of Path\newline' vario{v} ' variogram']);
    ylabel(['Standardized Frobenius Norm Error'])
end


figure(2);clf;hold on;
lim=8;
for v=1:numel(vario)
    for i=1:numel(nn)-1
        stat_eta{v}(i,:) = [mean(eta{v}{i}) std(eta{v}{i})];
        xx(i) = nn{i}(1);
        [f(i,:),xi(i,:)] = ksdensity(eta{v}{i});
        yi(i,:) = xx(i)*ones(1,numel(f(i,:)));
        ksdensity(eta{v}{i})
    end 
    %F=scatteredInterpolant(xi(:),yi(:),f(:))
    
    [X,Y]=meshgrid(linspace(min(min(xi)),max(max(xi)),100), linspace(min(min(yi)),max(max(yi)),256));
    
    
    a = F(X,Y);
    imagesc(X(1,:),Y(:,1),ones(256,100),'alphadata',a./max(max(a)))
   
    
%     plot(stat_eta{v}(:,1),xx)
%     plot(stat_eta{v}(:,1) - stat_eta{v}(:,2),xx,'--')
%     plot(stat_eta{v}(:,1) + stat_eta{v}(:,2),xx,'--')
end
ylabel('Number of Path');
xlabel('Standardized Frobenius Norm Error')
axis tight;

figure(1); clf; hold on;
lim=8;
filename2={'frobenium_sph_8_512', 'frobenium_sph_20_512' ,'frobenium_sph_52_512'};
col={'b','r','g','y'};
for i=1:numel(filename2)
    s=load(['result-SGSIM/Constant Path/' filename2{i}]);
    xlim_b=get(gca,'xlim');
    boxplot(cell2mat(s.eta(1:lim)),cell2mat(nn(1:lim)),'Orientation','horizontal','color',col{i},'symbol',['+' col{i}]);
    xlim_a=get(gca,'xlim');
    if i>1
        xlim([ min(xlim_b(1), xlim_a(1)) max(xlim_b(2), xlim_a(2))])
    end
end
hLegend = legend(findall(gca,'Tag','Box'), {'Neighboors: 8', 'Neighboors: 20', 'Neighboors: 52'});
ylabel('Number of Path');
xlabel('Standardized Frobenius Norm Error')


lim=8;
for v=1:numel(vario)
   (mean(eta{v}{5})-mean(eta{v}{1}))./mean(eta{v}{1})
end


