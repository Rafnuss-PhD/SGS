function [Sim] = definepath(Res,grid,parm)
Sim.xy=grid.xy(isnan(Res.m{1}));
Sim.n=numel(Sim.xy);

if parm.cstk
    Sim.path = definepath_in(Res,grid,parm,Sim);
    Sim.xy_r{1}=Sim.xy(Sim.path);
    [Sim.y_r{1},Sim.x_r{1}] = ind2sub([grid.ny, grid.nx],Sim.xy_r{1}); % * ind2sub() is taking linearized matrix index (i) and transform matrix index (i,j). We insert the randomized path of this specific simulation (ns) in Y.x and Y.y after the already known data of the primary data
else
    for i_realisation=1:parm.n_realisation
        Sim.path{i_realisation} = definepath_in(Res,grid,parm,Sim);
        Sim.xy_r{i_realisation}=Sim.xy(Sim.path{i_realisation});
        [Sim.y_r{i_realisation},Sim.x_r{i_realisation}] = ind2sub([grid.ny, grid.nx],Sim.xy_r{i_realisation}); % * ind2sub() is taking linearized matrix index (i) and transform matrix index (i,j). We insert the randomized path of this specific simulation (ns) in Y.x and Y.y after the already known data of the primary data
    end
end
end


function path = definepath_in(Res,grid,parm,Sim)
rng(parm.seed); % path among parrallel core can be the same if parm.seed is default
if strcmp(parm.path,'linear')
    if parm.path_random
        path = randperm(Sim.n);
    else
        path = 1:Sim.n;
    end
elseif  strcmp(parm.path,'spiralin') || strcmp(parm.path,'spiralout') || strcmp(parm.path,'maximize')
    dist_ = Inf*ones(Sim.n,1);
    X = Res.X(~isnan(Res.m{1}));
    Y = Res.Y(~isnan(Res.m{1}));
    for i_hard_data=1:numel(X)
        dist_chall = sqrt( (Res.X(Sim.xy)-X(i_hard_data)).^2 + (Res.Y(Sim.xy)-Y(i_hard_data)).^2 );
        dist_c = min(dist_,dist_chall(:));
    end
    if  strcmp(parm.path,'spiralout')
        if parm.path_random
            id_perm = randperm(numel(dist_));
            dist_perm = dist_(id_perm);
            [~,id_sort] = sort(dist_perm);
            path = id_perm(id_sort);
        else
            [~,path] = sort(dist_);
        end
    elseif strcmp(parm.path,'spiralin')
        if parm.path_random
            [dist_perm,id_perm] = randperm(dist_);
            [~,id_sort] = sort(dist_perm,'descend');
            path = id_perm(id_sort);
        else
            [~,path] = sort(dist_,'descend');
        end
    elseif strcmp(parm.path,'maximize')
        path = nan(Sim.n,1);
        for i_pt = 1:Sim.n
            if parm.path_random && i_pt~=1
                m=max(dist_);
                path(i_pt) = datasample(find(m==dist_),1); % randomize the selection in case of equal distence. (avoid a structural sampling)
            else
                [~,path(i_pt)] = max(dist_);
            end
            dist_chall = sqrt( (Res.X(Sim.xy)-Res.X(Sim.xy(path(i_pt)))).^2 + (Res.Y(Sim.xy)-Res.Y(Sim.xy(path(i_pt)))).^2 );
            dist_ = min(dist_,dist_chall);
        end
        
    end
    clear dist dist_chall X Y m
elseif strcmp(parm.path,'quasirandom')
    q = qrandstream('halton',2,'Skip',randi(1e10),'Leap',1e2);
    path =nan(Sim.n,1);
    i=1; j=1;
    while any(isnan(path))
        sub = floor(qrand(q,1).*[grid.nx grid.ny])+1;
        ind = sub2ind([grid.nx grid.ny], sub(1), sub(2) );
        if ~ismember(ind,Sim.xy(path(1:i-1))) && ismember(ind,Sim.xy)
            path(i) = find(Sim.xy==ind);
            i=i+1;
        end
        j=j+1;
    end
elseif numel(parm.path) == Sim.n
    parm.path = round(parm.path);
    assert(all(sort(parm.path(:)) == (1:Sim.n)'),'not valid path')
    path = round(parm.path);
else
    error('path method not define')
end

assert(numel(path)==numel(unique(path)),'not unique path')
end
