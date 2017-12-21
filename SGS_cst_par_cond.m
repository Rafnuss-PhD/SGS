%% Constant Path Sequential Gaussian Simulation
% See |SGS.m| for more general information on Sequential Gaussian
% Simulation. 
% SGS with a constant path uses a unique path for each realization, thus, it's code
% loop around the realization and the node after. See pseudo-code below
% <<fx1.jpg>>
% 
%%

function [Rest, t, parm] = SGS_cst_par_cond(nx,ny,m,covar,neigh,parm,hd)
tik.global = tic;

%% 1. Creation of the grid and path
[Y, X] = ndgrid(1:ny,1:nx);
tik.path = tic;
Path = nan(ny,nx);
Path(hd.id) = 0;
rng(parm.seed_path);
if parm.mg
   sx = 1:ceil(log(nx+1)/log(2));
   sy = 1:ceil(log(ny+1)/log(2));
   sn = max([numel(sy), numel(sx)]);
   nb = nan(sn,1);
   start = zeros(sn+1,1);
   ds = 2.^(sn-1:-1:0);
   path = nan(sum(isnan(Path(:))),1);
   for i_scale = 1:sn
       [Y_s,X_s] = ndgrid(1:ds(i_scale):ny,1:ds(i_scale):nx); % matrix coordinate
       id = find(isnan(Path(:)) & ismember([Y(:) X(:)], [Y_s(:) X_s(:)], 'rows'));
       nb(i_scale) = numel(id);
       start(i_scale+1) = start(i_scale)+nb(i_scale);
       path( start(i_scale)+(1:nb(i_scale)) ) = id(randperm(nb(i_scale)));
       Path(path( start(i_scale)+(1:nb(i_scale)) )) = start(i_scale)+(1:nb(i_scale));
       % Find the scaloe of hard data.
       hd.scale( ismember([hd.y hd.x], [Y_s(:) X_s(:)],'rows') & isnan(hd.scale)) =i_scale;
   end
else
   id=find(isnan(Path));
   path = id(randperm(numel(id)));
   Path(path) = 1:numel(id);
   ds=1; nb = numel(id); start=[0 nb]; sn=1;
end
t.path = toc(tik.path);


%% 2. Initialization Spiral Search
% Initialize spiral search stuff which don't change
x = ceil( min(covar(1).range(2)*neigh.wradius, nx));
y = ceil( min(covar(1).range(1)*neigh.wradius, ny));
[ss_Y, ss_X] = ndgrid(-y:y, -x:x);% grid{i_scale} of searching windows
ss_dist = sqrt( (ss_X/covar(1).range(2)).^2 + (ss_Y/covar(1).range(1)).^2); % find distence
ss_id_1 = find(ss_dist <= neigh.wradius); % filter node behind radius.
rng(parm.seed_search);
ss_id_1 = ss_id_1(randperm(numel(ss_id_1)));
[ss_dist_s, ss_id_2] = sort(ss_dist(ss_id_1)); % sort according distence.
ss_X_s=ss_X(ss_id_1(ss_id_2)); % sort the axis
ss_Y_s=ss_Y(ss_id_1(ss_id_2));
ss_n=numel(ss_X_s); %number of possible neigh

if parm.mg
    ss_scale=sn*ones(size(ss_X));
    for i_scale = sn-1:-1:1
        x_s = [-fliplr(ds(i_scale):ds(i_scale):x(end)) 0 ds(i_scale):ds(i_scale):x(end)]+(x+1);
        y_s = [-fliplr(ds(i_scale):ds(i_scale):y(end)) 0 ds(i_scale):ds(i_scale):y(end)]+(y+1);
        ss_scale(y_s,x_s)=i_scale;
    end
    ss_scale_s = ss_scale(ss_id_1(ss_id_2));
else
    ss_scale_s = sn*ones(size(ss_id_2));
end

%% 3. Initialization Covariance Lookup Table
if neigh.lookup
    ss_a0_C = zeros(ss_n,1);
    ss_ab_C = zeros(ss_n);
    for i=1:numel(covar)
        a0_h = sqrt(sum(([ss_Y_s(:) ss_X_s(:)]*covar(i).cx).^2,2));
        ab_h = squareform(pdist([ss_Y_s ss_X_s]*covar(i).cx));
        
        ss_a0_C = ss_a0_C + kron(covar(i).g(a0_h), covar(i).c0);
        ss_ab_C = ss_ab_C + kron(covar(i).g(ab_h), covar(i).c0);
    end
end


%% 4. Initizialization of the kriging weights and variance error
tik.weight = tic;
NEIGH = nan(nx*ny,neigh.nb);
% NEIGH_1 = nan(nx*ny,neigh.nb);
% NEIGH_2 = nan(nx*ny,neigh.nb);
LAMBDA = nan(nx*ny,neigh.nb);
S = nan(nx*ny,1);

XY_i=[Y(path) X(path)];

%% 5 Loop of scale for multi-grid path
for i_scale = 1:sn
    %% 5.1 Initializsed the search table of neighbors for the scale
    ss_id = find(ss_scale_s<=i_scale);
    ss_XY_s_s = [ss_Y_s(ss_id) ss_X_s(ss_id)];
    ss_dist_s_s = ss_dist_s(ss_id);
    if neigh.lookup
        ss_a0_C_s = ss_a0_C(ss_id);
        ss_ab_C_s = ss_ab_C(ss_id,ss_id);
    else
        ss_a0_C_s=[];
        ss_ab_C_s=[];
    end
    
    %% 5.2 Remove hard data which are on the current scale
    hd_XY_s = [hd.y(hd.scale>i_scale) hd.x(hd.scale>i_scale)];
    
    %% 5.3 Loop of simulated node
    for i_pt = start(i_scale)+(1:nb(i_scale))
        %% 5.3.1 Hard data
        hd_XY_d = bsxfun(@minus,hd_XY_s,XY_i(i_pt,:));
        hd_XY_d = hd_XY_d(hd_XY_d(:,1)<covar(1).range(1)*neigh.wradius &  hd_XY_d(:,2)<covar(1).range(2)*neigh.wradius,:);
        hd_dist=zeros(size(hd_XY_d,1),1);
        for i=1:numel(covar)
            hd_dist=hd_dist+sqrt(sum((hd_XY_d*covar(i).cx).^2,2));
        end
        
        [~, ss_hd_id] = sort( [ hd_dist; ss_dist_s_s]);
        tmp = [hd_XY_d; ss_XY_s_s];
        ss_hd_XY_s_s = tmp(ss_hd_id,:);
        
        %% 5.3.2 Neighborhood search
        n=0;
        neigh_nn=nan(neigh.nb,1);
        NEIGH_1 = nan(neigh.nb,1);
        NEIGH_2 = nan(neigh.nb,1);
        for nn = 2:size(ss_hd_XY_s_s,1) % 1 is the point itself... therefore unknown
            ijt = XY_i(i_pt,:) + ss_hd_XY_s_s(nn,:);
            if ijt(1)>0 && ijt(2)>0 && ijt(1)<=ny && ijt(2)<=nx
                if Path(ijt(1),ijt(2)) < i_pt % check if it,jt exist
                    n=n+1;
                    neigh_nn(n) = nn;
                    NEIGH_1(n) = ijt(1);
                    NEIGH_2(n) = ijt(2);
                    if n >= neigh.nb
                        break;
                    end
                end
            end
        end
        
        %% 5.3.3 Kriging system solving and storing of weights
        if n==0
            S(i_pt) = sum([covar.c0]);
        else
            NEIGH(i_pt,:) = NEIGH_1 + (NEIGH_2-1)* ny;
            if neigh.lookup
                a0_C = ss_a0_C_s(neigh_nn(1:n));
                ab_C = ss_ab_C_s(neigh_nn(1:n), neigh_nn(1:n));
            else
                D = pdist([0 0; ss_hd_XY_s_s(neigh_nn(1:n),:)]*covar.cx);
                C = covar.g(D);
                if n==1
                    a0_C = C;
                    ab_C = 1;
                else
                    a0_C = C(1:n)';
                    % Equivalent to : squareform(C(n+1:end));
                    ab_C = diag(ones(n,1))*0.5;
                    ab_C(tril(true(n),-1)) = C(n+1:end);
                    ab_C = ab_C + ab_C';
                end
            end
            l = ab_C \ a0_C;
            LAMBDA(i_pt,1:n) = l;
            S(i_pt) = sum([covar.c0]) - l'*a0_C;
        end
    end
    % disp(['scale: ' num2str(i_scale) '/' num2str(sn)])
end
t.weight = toc(tik.weight);

if parm.saveit
    filename=['result-SGS/SIM-', parm.name ,'_', datestr(now,'yyyy-mm-dd_HH-MM-SS'), '.mat'];
    mkdir('result-SGS/')
    save(filename, 'parm', 'nx','ny','start','nb', 'path', 'sn', 'k','NEIGH','S','LAMBDA')
end

%% 6. Realization loop
tik.real = tic;
Rest = nan(ny,nx,m);
parm_seed_U = parm.seed_U;
for i_real=1:m
    Res=nan(ny,nx);
    Res(hd.id)=hd.d;
    rng(parm_seed_U);
    U=randn(ny,nx);
    %% 6.1 Loop over scale and node for simulation
    for i_scale = 1:sn
        for i_pt = start(i_scale)+(1:nb(i_scale))
            n = ~isnan(NEIGH(i_pt,:));
            Res(path(i_pt)) = LAMBDA(i_pt,n)*Res(NEIGH(i_pt,n))' + U(i_pt)*sqrt(S(i_pt));
        end
    end
    Rest(:,:,i_real) = Res;
end

if parm.saveit
    filename=['result-SGS/SIM-', parm.name ,'_', datestr(now,'yyyy-mm-dd_HH-MM-SS'), '.mat'];
    mkdir('result-SGS/')
    save(filename, 'parm','nx','ny', 'Rest', 't', 'k','U')
end

t.real = toc(tik.real);
t.global = toc(tik.global);
end
