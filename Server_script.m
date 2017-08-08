clear all; close all
load('result-BSGS/GEN-Run_1_2017-05-07_14-37');
[Sigma.X,Sigma.Y]=meshgrid(Sigma.x, Sigma.y);
parm.k.covar = gen.covar;
parm.kernel=kern;

parm.unit='';
parm.nscore = 1;
parm.cstk = true;
parm.seed = 'shuffle';
parm.scale=[grid_gen.sx;grid_gen.sy]; % no multigrid
parm.saveit = false;
parm.k.method = 'sbss'; % sort, sbss (29.6) or smart minkmex
parm.k.quad = 0;
parm.k.wradius = 1.3;
parm.plot.kernel=0;
parm.plot.ns= 0;
parm.plot.krig=0;

parm.k.nb = [0 0 0 0 0; 10 10 10 10 20];

parm.n_realisation  = 0;

parm.par = 1;
parm.par_n = 4;

% use the log of hyd. cond.
Klog=K;
Klog.d=log(Klog.d);


% [Res, t, kern] = BSGS(Klog,Sigma,grid_gen,parm);

%% Weight

%% work on the weight
 filename = 'ResAB';
 parm.aggr.method='AB';
 parm.aggr.sum = 1;
 parm.aggr.a=0:1;%.01:.01:.19;
 parm.aggr.b=1/eps;
 [parm.aggr.A, parm.aggr.B] = meshgrid(parm.aggr.a, parm.aggr.b);
 parm.n_realisation  = parm.par_n*numel(parm.aggr.A);
% figure(929);hold on;
% i_pt_temp=1:grid_gen.nxy;
% fx = @(a,b,x) (atan(a*b) - atan(b*(a -  x )))/(atan(a*b) - atan(b*(a - 1)));
% for i_a=1:numel(parm.aggr.a)
%     a=parm.aggr.a(i_a);
%     for i_b=1:numel(parm.aggr.b)
%         b=parm.aggr.b(i_b);
%         plot(i_pt_temp,fx(a,b,i_pt_temp./grid_gen.nxy));
%     end
% end

%filename = 'ResSum155';
%parm.aggr.sum = 1.5;
%parm.aggr.method='cst';
%parm.n_realisation  = parm.par_n;
%parm.aggr.w=.5;


% filename = 'ResRAD';
% parm.aggr.method='rad';
% parm.aggr.a=[1 10 100];
% parm.aggr.b=[10 100 500 1000];
% [parm.aggr.A, parm.aggr.B] = meshgrid(parm.aggr.a, parm.aggr.b);
% parm.n_realisation  = parm.par_n*numel(parm.aggr.A);

% rad=Sigma.rad(randperm(length(Sigma.rad)));
% figure; hold on;
% for i_a=1:numel(parm.aggr.a)
%     a=parm.aggr.a(i_a);
%     for i_b=1:numel(parm.aggr.b)
%         b=parm.aggr.b(i_b);
%         w = 1-rad.*(a-b.*linspace(0,1,numel(rad)));
%         w(w<0) = 0;
%         w(w>1) = 1;
%         plot(w)
%     end
% end

 
disp('Setup ok, Start')
[Res,~,kern,~,~,Nscore] =  BSGS(Klog,Sigma,grid_gen,parm);
disp('Run finished')

% Sigma_d = Sigma.d(any( bsxfun(@eq,Sigma.X(:),Res.X(:)') & bsxfun(@eq,Sigma.Y(:),Res.Y(:)' ),2));
Gamma_t = (1-parm.k.covar(1).g(grid_gen.x/parm.k.covar(1).range(1)))';
id = Res.x<parm.k.covar(1).range0(1).*parm.k.wradius;
Gamma_t_id = Gamma_t( any( bsxfun(@eq, grid_gen.x,Res.x(id)) ));

E1 = nan(sum(id),parm.n_realisation);
E2 = nan(numel(parm.kernel.dens),parm.n_realisation);
    
XY = kern.XY;
Sigma_d = Sigma.d(:);
Res_m_ns = Res.m_ns;
Res_m = Res.m;
dens = parm.kernel.dens(:);
parfor i_realisation=1:parm.n_realisation
   gamma_x = variogram_gridded_perso(Res_m_ns{i_realisation});
   E1(:,i_realisation) = gamma_x(id)-Gamma_t_id;
   E2(:,i_realisation) = ksdensity([Sigma_d Res_m{i_realisation}(:)],XY)-dens;
end

try 
    if strcmp(parm.aggr.method,'AB')
        parfor i_ab = 1:numel(parm.aggr.A)
            idAB=i_ab:numel(parm.aggr.A):parm.n_realisation;
            OF1(i_ab) = sqrt(mean(mean(E1(:,idAB),2).^2));
            OF2(i_ab) = sqrt(mean(mean(E2(:,idAB),2).^2));
        end
    elseif strcmp(parm.aggr.method,'rad')
        parfor i_x = 1:numel(parm.aggr.x)
            idx=i_x:numel(parm.aggr.x):parm.n_realisation;
            OF1(i_x) = sqrt(mean(mean(E1(:,idx),2).^2));
            OF2(i_x) = sqrt(mean(mean(E2(:,idx),2).^2));
        end
    else
        OF1 = sqrt(mean(mean(E1,2).^2));
        OF2 = sqrt(mean(mean(E2,2).^2));
    end
catch
    disp('Error for computing OF')
    save(['./result-BSGS/' filename],'Res','parm','kern','Nscore','E1','E2')
end

save(['./result-BSGS/' filename],'Res','parm','kern','Nscore','E1','E2','OF1','OF2')
disp('file written')
%save('result-BSS/ResAB_A_1.mat','parm','E1','E2')
% save('C:\Users\Raphael\Desktop\to be compress\ResRAD_AB.mat','Res','parm','kern','Nscore','E1','E2','OF1','OF2')


