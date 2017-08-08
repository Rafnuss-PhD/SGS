function [pt_krig] = neighsearch(pt,Res,Prim,krig,parm,i_realisation)

if strcmp(krig.method,'smart')
    warning('Smart Search not good yet!')
   if pt.i<Res.sim.n/3
        krig.method='sort';
   else
       krig.method='sbss';
   end
end


if strcmp(krig.method,'sbss')
    % 1. Super Grid Block from Hard Data:
    sb_j = min([round((Res.y(pt.y)-krig.sb.y(1))/krig.sb.dy +1)'; krig.sb.ny]);
    sb_i = min([round((Res.x(pt.x) -krig.sb.x(1))/krig.sb.dx +1)'; krig.sb.nx]);
    mask.prim_all = reshape(krig.sb.mask(sb_j,sb_i,:),size(krig.sb.mask,3),1);
    
    [~,sb_id] = sort( sqrt( ( (Prim.x(mask.prim_all)-Res.x(pt.x))./krig.covar(1).range(1) ).^2 + ( (Prim.y(mask.prim_all)-Res.y(pt.y))./krig.covar(1).range(2) ).^2 ) );
    mask.prim_all_id = find(mask.prim_all);
    pt_krig.mask.prim = false(Prim.n,1);
    pt_krig.mask.prim( mask.prim_all_id(sb_id(1:min(krig.nb(2),numel(sb_id)))))=true;
    
    % 2. Spiral search
    if krig.quad % 2. per quandrant
        nn_max=length(krig.ss.el.dist_s);
        n=[0 0 0 0];
        sel_ss_idx=cell(4,1);
        for q=1:4
            sel_ss_idx{q}=nan(krig.nb(2,q),2);
            nn=2; % 1 is the point itself... therefore unknown
            while n(q)<krig.nb(2,q) && nn<=nn_max && krig.ss.el.dist_s(nn)<=krig.wradius % while not exceed number of point wanted and still inside the ellipse
                it = pt.x + krig.qs(q,1)*krig.ss.el.X_s(nn);
                jt = pt.y + krig.qs(q,2)*krig.ss.el.Y_s(nn);
                if it>0 && jt>0 && it<=Res.nx && jt <=Res.ny  && all(Prim.x(pt_krig.mask.prim)~=it) && all(Prim.y(pt_krig.mask.prim)~=jt) % check to not be outside the grid
                    if ~isnan(Res.m_ns{i_realisation}(jt,it)) % check if it,jt exist
                        n(q)=n(q)+1;
                        sel_ss_idx{q}(n(q),:) = [jt it];
                    end
                end
                nn=nn+1;
            end
            sel_ss_idx{q}=sel_ss_idx{q}(1:n(q),:); % only the max number of point found.
        end
        k0_ss_idx = unique([sel_ss_idx{1};sel_ss_idx{2};sel_ss_idx{3};sel_ss_idx{4}],'rows');
        pt_krig.mask.res = sub2ind([Res.ny, Res.nx],k0_ss_idx(:,1),k0_ss_idx(:,2));
    else
        n=0;
        sel_ss_idx=nan(krig.n,2);
        nn=2; % 1 is the point itself... therefore unknown

        while n<krig.n && nn<=krig.ss.el.n % while not exceed number of point wanted and still inside the ellipse
            it = pt.x + krig.ss.el.X_s(nn);
            jt = pt.y + krig.ss.el.Y_s(nn);
            if it>0 && jt>0 && it<=Res.nx && jt <=Res.ny %&&  ~any(Prim.x(pt_krig.mask.prim)==Res.x(it) & Prim.y(pt_krig.mask.prim)==Res.y(jt)) % check to not be outside the grid
                if ~isnan(Res.m_ns{i_realisation}(jt,it)) % check if it,jt exist
                    n=n+1;
                    sel_ss_idx(n,:) = [jt it];
                end
            end
            nn=nn+1;
        end
        sel_ss_idx=sel_ss_idx(1:n,:); % only the max number of point found.
    end
    pt_krig.mask.res = sub2ind([Res.ny, Res.nx],sel_ss_idx(:,1),sel_ss_idx(:,2));

    
elseif strcmp(krig.method,'sort')
    mask.res = Res.sim.xy_r{i_realisation}(1:pt.i-1);
    sel_g_ini=[ [Prim.x Prim.y]; [Res.x(Res.sim.x_r{i_realisation}(1:pt.i-1))  Res.y(Res.sim.y_r{i_realisation}(1:pt.i-1))]];
    % Just remove point outside search radius and keep
    % This is identical (more or less) to cokri_cam (sort and selection)
    center = [Res.x(pt.x) Res.y(pt.y)];
    dist = sqrt( ((sel_g_ini(:,1)-center(1))/krig.covar(1).range(1)).^2 + ((sel_g_ini(:,2)-center(2))/krig.covar(1).range(2)).^2 );
    [dist_s, idx_s] = sort(dist);
    
    if krig.quad
        % Compute angle of all neigh from sim pt: define their beloging to
        % quadrant using the range of radius
        dist_angl = atan2(sel_g_ini(:,1)-center(1), sel_g_ini(:,2)-center(2));
        dist_angl = dist_angl(idx_s);
        dist_cat = false(numel(dist_s),4);
        dist_cat(0<=dist_angl&dist_angl<=pi/2,1) = true;
        dist_cat(pi/2<=dist_angl&dist_angl<=pi,2) = true;
        dist_cat(-pi<=dist_angl&dist_angl<=-pi/2,3) = true;
        dist_cat(-pi/2<=dist_angl&dist_angl<=0,4) = true;
        sb_i=1;
        quad=randperm(4); % permutate the order of quadrant to avoid preferentially one side
        quad_i=1;
        mask.all=[];
        
        while sb_i<=sum(krig.nb(2,:)) && numel(dist_s)>0 && any(dist_s < krig.wradius)
            sb_id = find(dist_cat(:,quad(quad_i)), 1, 'first'); % find the closest point in the quandrant
            if ~isempty(sb_id) % if there is a point in the quandrant, put in in the list to keep
                mask.all=[mask.all; idx_s(sb_id)];
                dist_cat(sb_id,:)=[];
                dist_s(sb_id)=[];
                idx_s(sb_id)=[];
                sb_i=sb_i+1;
            else % otherwise, delete the quandrant in order to not to look anymore inside this quadrant
                quad(quad_i)=[];
                quad_i=quad_i-1;
            end
            quad_i = max(1,mod(quad_i+1,numel(quad)+1)); % next wuadrant (mod(...) return to 1 once looped.
        end
        
    else
        sb_i=1;
        mask.all=nan(sum(krig.nb(2,:)),1);
        while sb_i<=sum(krig.nb(2,:)) && sb_i<=numel(dist_s) && dist_s(sb_i) < krig.wradius
            mask.all(sb_i)=idx_s(sb_i);
            sb_i=sb_i+1;
        end
        mask.all = mask.all(~isnan(mask.all));
        % if numel(mask.all)==0;
        %     mask.all=zeros(0,1);
        % end
    end
    
    pt_krig.mask.prim = mask.all(mask.all<=Prim.n);
    pt_krig.mask.res = mask.res(mask.all(mask.all>Prim.n)-Prim.n);
    
    % Force hard data
    %     if ~ismember([32 32], sel_g,'rows')
    %         sel_g=[ sel_g(1:end-1,:); 32 32];
    %         pt_krig.mask = [pt_krig.mask(1:end-1); find(ismember(sel_g_ini,[32 32],'rows')) ];
    %     end
    
    

elseif strcmp(krig.method,'minkmex')
    
    mask.res = Res.sim.xy_r{i_realisation}(1:pt.i-1);
    sel_g_ini=[ [Prim.x Prim.y]; [Res.x(Res.sim.x_r{i_realisation}(1:pt.i-1))  Res.y(Res.sim.y_r{i_realisation}(1:pt.i-1))]];

    center = [Res.x(pt.x) Res.y(pt.y)];
    dist = sqrt( ((sel_g_ini(:,1)-center(1))/krig.covar(1).range(1)).^2 + ((sel_g_ini(:,2)-center(2))/krig.covar(1).range(2)).^2 );

    mask.all = find(dist <= krig.wradius);
    % idx = partial_sort_index(dist(mask.all), sum(krig.nb(2,:)));
    [~, idx] = minkmex(dist(mask.all), sum(krig.nb(2,:)));
    mask.all = mask.all(idx);
    
    pt_krig.mask.prim = mask.all(mask.all<=Prim.n);
    pt_krig.mask.res = mask.res(mask.all(mask.all>Prim.n)-Prim.n);
    

elseif strcmp(krig.method,'partialSort')
   
    mask.res = Res.sim.xy_r{i_realisation}(1:pt.i-1);
    sel_g_ini=[ [Prim.x Prim.y]; [Res.x(Res.sim.x_r{i_realisation}(1:pt.i-1))  Res.y(Res.sim.y_r{i_realisation}(1:pt.i-1))]];

    center = [Res.x(pt.x) Res.y(pt.y)];
    dist = sqrt( ((sel_g_ini(:,1)-center(1))/krig.covar(1).range(1)).^2 + ((sel_g_ini(:,2)-center(2))/krig.covar(1).range(2)).^2 );

    mask.all = find(dist <= krig.wradius);
    % idx = partial_sort_index(dist(mask.all), sum(krig.nb(2,:)));
    if numel(mask.all)>0
        [~, idx] = firstKSmallest(dist(mask.all), sum(krig.nb(2,:)));
        idx=idx';
    else
        idx=mask.all;
    end

    mask.all = mask.all(idx);
    
    pt_krig.mask.prim = mask.all(mask.all<=Prim.n);
    pt_krig.mask.res = mask.res(mask.all(mask.all>Prim.n)-Prim.n);
    
elseif strcmp(krig.method,'knnsearch')
    mask.res = Res.sim.xy_r{i_realisation}(1:pt.i-1);
    sel_g_ini=[ [Prim.x Prim.y]; [Res.x(Res.sim.x_r{i_realisation}(1:pt.i-1))  Res.y(Res.sim.y_r{i_realisation}(1:pt.i-1))]];

    center = [Res.x(pt.x) Res.y(pt.y)];
   
    IDX = knnsearch(sel_g_ini,center,'K',krig.n);
    
    mask.all = IDX;    
    
    pt_krig.mask.prim = mask.all(mask.all<=Prim.n);
    pt_krig.mask.res = mask.res(mask.all(mask.all>Prim.n)-Prim.n);
        
else
    error(['Not possible neigh search method: ' krig.method ' (smart,sbss,nearest)'])
end
