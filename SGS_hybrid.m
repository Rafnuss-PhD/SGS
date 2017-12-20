function [Rest, t, k, parm] = SGS_hybrid(ny,nx,parm)

tik.global = tic;


%% 1. Creation of the grid and path
tik.path = tic;
[Y, X] = ndgrid(1:ny,1:nx);
Path = nan(ny,nx);
rng(parm.seed_path);
sx = 1:ceil(log(nx+1)/log(2));
sy = 1:ceil(log(ny+1)/log(2));
sn = max([numel(sy), numel(sx)]);
ds = 2.^(sn-1:-1:0);
nb = nan(sn,1);
start = zeros(sn+1,1);
path = nan(nx*ny,1);
if (parm.seed_path>1)
    [Y_s,X_s] = ndgrid(1:ds(parm.seed_path-1):ny,1:ds(parm.seed_path-1):nx); % matrix coordinate
    Path(ismember([Y(:) X(:)], [Y_s(:) X_s(:)], 'rows')) = 0;
    nb(parm.seed_path-1) = sum(Path(:)==0);
    start(parm.seed_path) = nb(parm.seed_path-1);
end
for i_scale = parm.seed_path:sn
    [Y_s,X_s] = ndgrid(1:ds(i_scale):ny,1:ds(i_scale):nx); % matrix coordinate
    id = find( isnan(Path(:)) & ismember([Y(:) X(:)], [Y_s(:) X_s(:)], 'rows'));
    nb(i_scale) = numel(id);
    start(i_scale+1) = start(i_scale)+nb(i_scale);
    path( start(i_scale)+(1:nb(i_scale)) ) = id(randperm(nb(i_scale)));
    Path(path( start(i_scale)+(1:nb(i_scale)) )) = start(i_scale)+(1:nb(i_scale));
end
t.path = toc(tik.path);


%% 2. Initialization Spiral Search
% Initialize spiral search stuff which don't change
x = ceil( min(k.covar(1).range(2)*k.wradius, nx));
y = ceil( min(k.covar(1).range(1)*k.wradius, ny));
[ss_Y, ss_X] = ndgrid(-y:y, -x:x);% grid{i_scale} of searching windows
ss_dist = sqrt( (ss_X/k.covar(1).range(2)).^2 + (ss_Y/k.covar(1).range(1)).^2); % find distence
ss_id_1 = find(ss_dist <= k.wradius); % filter node behind radius.
rng(parm.seed_search);
ss_id_1 = ss_id_1(randperm(numel(ss_id_1)));
[~, ss_id_2] = sort(ss_dist(ss_id_1)); % sort according distence.
ss_X_s=ss_X(ss_id_1(ss_id_2)); % sort the axis
ss_Y_s=ss_Y(ss_id_1(ss_id_2));
ss_n=numel(ss_X_s); %number of possible neigh

ss_scale=sn*ones(size(ss_X));
for i_scale = sn-1:-1:1
    x_s = [-fliplr(ds(i_scale):ds(i_scale):x(end)) 0 ds(i_scale):ds(i_scale):x(end)]+(x+1);
    y_s = [-fliplr(ds(i_scale):ds(i_scale):y(end)) 0 ds(i_scale):ds(i_scale):y(end)]+(y+1);
    ss_scale(y_s,x_s)=i_scale;
end
ss_scale_s = ss_scale(ss_id_1(ss_id_2));


%% 3. Initialization of covariance lookup table
tik.covtable = tic;
ss_a0_C = zeros(ss_n,1);
ss_ab_C = zeros(ss_n);
for i=1:numel(k.covar)
    a0_h = sqrt(sum(([ss_Y_s(:) ss_X_s(:)]*k.covar(i).cx).^2,2));
    ab_h = squareform(pdist([ss_Y_s ss_X_s]*k.covar(i).cx));
    ss_a0_C = ss_a0_C + kron(k.covar(i).g(a0_h), k.covar(i).c0);
    ss_ab_C = ss_ab_C + kron(k.covar(i).g(ab_h), k.covar(i).c0);
end
% Transform ss.ab_C sparse?
k_nb = neigh.nb;
k_covar_c0 = sum([k.covar.c0]);
t.covtable = toc(tik.covtable);

%% Constant path Loop
tik.weight = tic;
NEIGH = nan(nx*ny,neigh.nb);
% NEIGH_1 = nan(nx*ny,neigh.nb);
% NEIGH_2 = nan(nx*ny,neigh.nb);
LAMBDA = nan(nx*ny,neigh.nb);
S = nan(nx*ny,1);

XY_i=nan(nx*ny,2);
XY_i(~isnan(path),:) = [Y(path(~isnan(path))) X(path(~isnan(path)))];


for i_scale = parm.seed_path:sn
    ss_id = find(ss_scale_s<=i_scale);
    ss_XY_s = [ss_Y_s(ss_id) ss_X_s(ss_id)];
    ss_a0_C_s = ss_a0_C(ss_id);
    ss_ab_C_s = ss_ab_C(ss_id,ss_id);
    
    for i_pt = start(i_scale)+(1:nb(i_scale))
        n=0;
        neigh_nn=nan(k_nb,1);
        NEIGH_1 = nan(k_nb,1);
        NEIGH_2 = nan(k_nb,1);
        for nn = 2:size(ss_XY_s,1) % 1 is the point itself... therefore unknown
            ijt = XY_i(i_pt,:) + ss_XY_s(nn,:);
            if ijt(1)>0 && ijt(2)>0 && ijt(1)<=ny && ijt(2)<=nx
                if Path(ijt(1),ijt(2)) < i_pt % check if it,jt exist
                    n=n+1;
                    neigh_nn(n) = nn;
                    NEIGH_1(n) = ijt(1);
                    NEIGH_2(n) = ijt(2);
                    if n >= k_nb
                        break;
                    end
                end
            end
        end
        if n==0
            S(i_pt) = k_covar_c0;
        else
            NEIGH(i_pt,:) = NEIGH_1 + (NEIGH_2-1)* ny;
             a0_C = ss_a0_C_s(neigh_nn(1:n));
            ab_C = ss_ab_C_s(neigh_nn(1:n), neigh_nn(1:n));
            l = ab_C \ a0_C;
            LAMBDA(i_pt,1:n) = l;
            S(i_pt) = k_covar_c0 - l'*a0_C;
        end
    end
end
t.weight = toc(tik.weight);


%% Realization loop
tik.real = tic;
Rest = nan(ny,nx,parm.n_real);

for i_real=1:parm.n_real
    Res=nan(ny,nx);
    % Path
    Path_real = Path;
    path_real = path;
    rng(parm.seed_path);
    for i_scale = 1:parm.seed_path-1
        [Y_s,X_s] = ndgrid(1:ds(i_scale):ny,1:ds(i_scale):nx); % matrix coordinate
        id = find(Path_real(:)==0 & ismember([Y(:) X(:)], [Y_s(:) X_s(:)], 'rows'));
        nb(i_scale) = numel(id);
        start(i_scale+1) = start(i_scale)+nb(i_scale);
        path_real( start(i_scale)+(1:nb(i_scale)) ) = id(randperm(nb(i_scale)));
        Path_real(path_real( start(i_scale)+(1:nb(i_scale)) )) = start(i_scale)+(1:nb(i_scale));
    end
    XY_i=[Y(path_real) X(path_real)];
    rng(parm.seed_U);
    U=randn(ny,nx);
    for i_scale = 1:parm.seed_path-1
        % Neighborhood Search
        ss_id = find(ss_scale_s<=i_scale);
        ss_XY_s = [ss_Y_s(ss_id) ss_X_s(ss_id)];
        ss_a0_C_s = ss_a0_C(ss_id);
        ss_ab_C_s = ss_ab_C(ss_id,ss_id);
        for i_pt = start(i_scale)+(1:nb(i_scale))
            n=0;
            neigh_nn=nan(k_nb,1);
            NEIGH_1 = nan(k_nb,1);
            NEIGH_2 = nan(k_nb,1);
            for nn = 2:size(ss_XY_s,1) % 1 is the point itself... therefore unknown
                ijt = XY_i(i_pt,:) + ss_XY_s(nn,:);
                if ijt(1)>0 && ijt(2)>0 && ijt(1)<=ny && ijt(2)<=nx
                    if Path_real(ijt(1),ijt(2)) < i_pt % check if it,jt exist
                        n=n+1;
                        neigh_nn(n) = nn;
                        NEIGH_1(n) = ijt(1);
                        NEIGH_2(n) = ijt(2);
                        if n >= k_nb
                            break;
                        end
                    end
                end
            end
            
            if n==0
                Res(path_real(i_pt)) = U(i_pt)*k_covar_c0;
            else
                a0_C = ss_a0_C_s(neigh_nn(1:n));
                ab_C = ss_ab_C_s(neigh_nn(1:n), neigh_nn(1:n));
                LAMBDA_t = ab_C \ a0_C;
                S_t = k_covar_c0 - LAMBDA_t'*a0_C;
                NEIGH_t = NEIGH_1 + (NEIGH_2-1)* ny;
                Res(path_real(i_pt)) = LAMBDA_t'*Res(NEIGH_t(1:n)) + U(i_pt)*sqrt(S_t);
            end
        end
        for i_scale = parm.seed_path:sn
            for i_pt = start(i_scale)+(1:nb(i_scale))
                n = ~isnan(NEIGH(i_pt,:));
                Res(path_real(i_pt)) = LAMBDA(i_pt,n)*Res(NEIGH(i_pt,n))' + U(i_pt)*sqrt(S(i_pt));
            end
        end
    end
    Rest(:,:,i_real) = Res;
end

t.real = toc(tik.real);
t.global = toc(tik.global);
end
