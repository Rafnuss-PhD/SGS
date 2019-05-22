%% Variance Covariance Matrix of Sequential Gaussian Simulation
% See |SGS.m| for more general information on Sequential Gaussian
% Simulation. 
% 
%%


function CY = SGS_varcovar(nx,ny,~,covar,neigh,parm)
%% 1. Creation of the grid and path
[Y, X] = ndgrid(1:ny,1:nx);
Path = nan(ny,nx);
rng(parm.seed_path);
if parm.mg
   sx = 1:ceil(log(nx+1)/log(2));
   sy = 1:ceil(log(ny+1)/log(2));
   sn = max([numel(sy), numel(sx)]);
   nb = nan(sn,1);
   start = zeros(sn+1,1);
   ds = 2.^(sn-1:-1:0);
   path = nan(nx*ny,1);
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


%% 2. Initialization Spiral Search
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
    ss_scale_s = sn*ones(size(ss_id_2));
end

%% 3. Initialization Covariance Lookup Table
ss_a0_C = zeros(ss_n,1);
ss_ab_C = zeros(ss_n);
for i=1:numel(covar)
    a0_h = sqrt(sum(([ss_Y_s(:) ss_X_s(:)]*covar(i).cx).^2,2));
    ab_h = squareform(pdist([ss_Y_s ss_X_s]*covar(i).cx));
    
    ss_a0_C = ss_a0_C + kron(covar(i).g(a0_h), covar(i).c0);
    ss_ab_C = ss_ab_C + kron(covar(i).g(ab_h), covar(i).c0);
end

%% 4. Initialization of variable               
LambdaM = zeros(nx*ny,neigh.nb+1);
XY_i=[Y(path) X(path)];

%% 5 Loop of scale for multi-grid path
for i_scale = 1:sn
    %% 5.1 Initializsed the search table of neighbors for the scale
    ss_id = find(ss_scale_s<=i_scale);
    ss_XY_s = [ss_Y_s(ss_id) ss_X_s(ss_id)];
    ss_a0_C_s = ss_a0_C(ss_id);
    ss_ab_C_s = ss_ab_C(ss_id,ss_id);
    
    %% 5.2 Loop of simulated node
    for i_pt = start(i_scale)+(1:nb(i_scale))
        %% 5.2.1 Neighborhood search
        n=0;
        neigh_nn=nan(neigh.nb,3);
        for nn = 2:size(ss_XY_s,1) % 1 is the point itself... therefore unknown
            ijt = XY_i(i_pt,:)+ss_XY_s(nn,:);
            if ijt(1)>0 && ijt(2)>0 && ijt(1)<=ny && ijt(2)<=nx 
                if Path(ijt(1),ijt(2)) < i_pt % check if it,jt exist
                    n=n+1;
                    neigh_nn(n,:) = [nn ijt];
                    if n >= neigh.nb
                        break;
                    end
                end
            end
        end
        
        %% 5.2.2 Kriging system solving and storing of weights
        neigh_id = (neigh_nn(1:n,2:3)-1)*[1 ny]'+1;
        a0_C = ss_a0_C_s(neigh_nn(1:n,1));
        ab_C = ss_ab_C_s(neigh_nn(1:n,1), neigh_nn(1:n,1));
        l = ab_C \ a0_C;
        S = sum([covar.c0]) - l'*a0_C;
        LambdaM(path(i_pt), neigh_id) = l./sqrt(S);
        LambdaM(path(i_pt), path(i_pt)) = -1/sqrt(S);
    end
end

%% 6 Compute the Covariance 
CY = sparse(LambdaM) \ transpose(inv(sparse(LambdaM)));
end
