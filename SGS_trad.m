%% Traditional Sequential Gaussian Simulation
% See |SGS.m| for more general information on Sequential Gaussian
% Simulation. 
% Traditional SGS uses a new path for each realization, thus, it's code
% loop around the realization and the node after. See pseudo-code below
% <<fx1.jpg>>
% 
%%

function [Rest, t] = SGS_trad(nx,ny,m,covar,neigh,parm)
tik.global = tic;

%% 1. Creation of the grid and initialization of path
[Y, X] = ndgrid(1:ny,1:nx);
if parm.mg
    sx = 1:ceil(log(nx+1)/log(2));
    sy = 1:ceil(log(ny+1)/log(2));
    sn = max([sy(end), sx(end)]);
    nb = nan(sn,1);
    start = zeros(sn+1,1);
    ds = 2.^(sn-1:-1:0);
end


%% 2. Initialization Spiral Search
% Initialize spiral search stuff which don't change
x = ceil( min(covar(1).range(2)*neigh.wradius, nx));
y = ceil( min(covar(1).range(1)*neigh.wradius, ny));
[ss_Y, ss_X] = ndgrid(-y:y, -x:x);% grid{i_scale} of searching windows
ss_dist = sqrt( (ss_X/covar(1).range(2)).^2 + (ss_Y/covar(1).range(1)).^2); % find distence
ss_id_1 = find(ss_dist <= neigh.wradius); % filter node behind radius.
rng(parm.seed_search);
ss_id_1 = ss_id_1(randperm(numel(ss_id_1)));
[~, ss_id_2] = sort(ss_dist(ss_id_1)); % sort according distence.
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
    ss_scale_s = ones(size(ss_id_2));
end


%% 3. Initialization of covariance lookup table
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


%% 4. Realization loop
Rest = nan(ny,nx,m);
tik.real = tic;

for i_real=1:m
    Res=nan(ny,nx);
    
    %% 4.1 Generation of the path
    Path = nan(ny,nx);
    path = nan(nx*ny,1);
    rng(parm.seed_path);
    if parm.mg
        for i_scale = 1:sn
            [Y_s,X_s] = ndgrid(1:ds(i_scale):ny,1:ds(i_scale):nx); % matrix coordinate
            id = find(isnan(Path(:)) & ismember([Y(:) X(:)], [Y_s(:) X_s(:)], 'rows'));
            nb(i_scale) = numel(id);
            start(i_scale+1) = start(i_scale)+nb(i_scale);
            path( start(i_scale)+(1:nb(i_scale)) ) = id(randperm(nb(i_scale)));
            Path(path( start(i_scale)+(1:nb(i_scale)) )) = start(i_scale)+(1:nb(i_scale));
        end
    else
        id=find(isnan(Path));
        path = id(randperm(numel(id)));
        Path(path) = 1:numel(id);
        ds=1; nb = numel(id); start=[0 nb]; sn=1;
    end
    XY_i=[Y(path) X(path)];
    
    rng(parm.seed_U);
    U=randn(ny,nx);
    
    %% 4.2 Loop of scale for multi-grid path
    for i_scale = 1:sn
        
        %% 4.2.1 Initializsed the search table of neighbors for the scale
        ss_id = find(ss_scale_s<=i_scale);
        ss_XY_s = [ss_Y_s(ss_id) ss_X_s(ss_id)];
        if neigh.lookup
            ss_a0_C_s = ss_a0_C(ss_id);
            ss_ab_C_s = ss_ab_C(ss_id,ss_id);
        end
        
        %% 4.2.1 Loop of simulated node
        for i_pt = start(i_scale)+(1:nb(i_scale))
            %% 4.2.1.1 Neighborhood search
            n=0;
            neigh_nn=nan(neigh.nb,1);
            NEIGH_1 = nan(neigh.nb,1);
            NEIGH_2 = nan(neigh.nb,1);
            for nn = 2:size(ss_XY_s,1) % 1 is the point itself... therefore unknown
                ijt = XY_i(i_pt,:) + ss_XY_s(nn,:);
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
            
            %% 4.2.1.2 Kriging system solving and simulation
            if n==0
                Res(path(i_pt)) = U(i_pt)*sum([covar.c0]);
                NEIGH=[];
            else
                if neigh.lookup
                    a0_C = ss_a0_C_s(neigh_nn(1:n));
                    ab_C = ss_ab_C_s(neigh_nn(1:n), neigh_nn(1:n));
                else
                    D = pdist([0 0; ss_XY_s(neigh_nn(1:n),:)]*covar.cx);
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
                
                LAMBDA = ab_C \ a0_C;
                S = sum([covar.c0]) - LAMBDA'*a0_C;
                
                NEIGH = NEIGH_1 + (NEIGH_2-1)* ny;
                Res(path(i_pt)) = LAMBDA'*Res(NEIGH(1:n)) + U(i_pt)*sqrt(S);
            end
        end
    end
    Rest(:,:,i_real) = Res;
end

t.real = toc(tik.real);
t.global = toc(tik.global);
end
