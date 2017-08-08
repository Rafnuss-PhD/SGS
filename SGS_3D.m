%% Sequential Gaussian Simulation (SGS)
% tags: 2D, Constant Path, Covariance Matrix, parralele computing
%
%
% INPUT:
%
% * Prim.d      : Primary variable (1D-X.n elts)
% * Prim.x,y    : x,y-coordinate of the primary variable (1D-X.n elts)
% * parm        : Parameters of the simulation.
%   * scale     : array of the scales (i, where nx=2^grid{i}+1) to simulate
%   * seed      : random number. use to reproduce exactly the same simulation.
%   * saveit 	: save or not in a file
%   * name 	    : filename to be save
%   * familyname: folder where to save the file
%   * unit 	    : unit of the primary variable, used for plot
%   * n_realisation : number of realisations
%   * par       : bolean. parralize relization
%   * par_n     : specific the number of core
%   * dist      : 
%   * notify    : bolean. send email
%   * notify_email : string email
%   * nscore	: bolean. with or without  normal score transform
%   * path      : type of path
%   * path_random : bolean. with or without randomized the seed of path
%   * varcovar  : 
%   * cstk      : bolean. with or without constant path
%   * cstk_s 	: scale at which cst path is switch on
%   * k:
%       * nb_neigh: row: min and max number of point to use in each
%       quadrant (first 4 columns) and hard data (5th col.)
%       * method: 
%       * sb:
%           * nx:
%           * ny:
%       * wradius (1.3):
%   * plot:
%       * bsgs	: plot. ?
%       * ns 	: plot. ?
%       * sb 	: plot. ?
%       * kernel : plot. ?
%       * fitvar : plot. ?
%       * krig  : plot. ?
%
%
% OUTPUT:
%
% * Res     : Resultant variable
% * t       : Time of simulations
% * k       : kriging information
% * parm    : Parameter structure
% * filename: Name of the file save
%
%
% * *Author:* Raphael Nussbaumer (raphael.nussbaumer@unil.ch)
%

function [Res, t, k, grid.g, parm, filename] = SGS(Z,parm)
t.global = tic;
% addpath(genpath('./functions/.'))

%% * *INPUT CEHCKING*
% This section of the code generates a valid parm structure based on the
% inputted parm. If some parm value haven't been input, this section will
% fill automatically with defautl value. This may not allowed be the best.



% Parameter settings
if ~isfield(parm, 'seed'),          parm.seed           = rand(); end
rng(parm.seed);
if ~isfield(parm, 'saveit'),        parm.saveit         = 1; end % bolean, save or not the result of simulation
if ~isfield(parm, 'name'),          parm.name           = ''; end % name use for saving file
if ~isfield(parm, 'familyname'),    parm.familyname     = ''; end
if ~isfield(parm, 'unit'),          parm.unit           = ''; end % unit of the primary variable, used in plot
if ~isfield(parm, 'n_realisation'), parm.n_realisation  = 1; end
if ~isfield(parm, 'par'),           parm.par            = 1; end
if ~isfield(parm, 'dist'),          parm.dist           = 0; end
if ~isfield(parm, 'par_n'),         parm.par_n          = feature('numcores'); end
if ~isfield(parm, 'notify'),
    parm.notify          = 0;
else
    if ~isfield(parm, 'notify_email'), parm.notify_email  = 'rafnuss@gmail.com'; end
end

% Path
if ~isfield(parm, 'path'),          parm.path            = 'linear'; end
if ~isfield(parm, 'path_random'),   parm.path_random     = true; end
if ~isfield(parm, 'varcovar'),      parm.path            = 0; end


% Grid
if ~isfield(grid, 'scale')
    grid.scale = [6; 6];
end
[grid.dim, grid.n_scale] = size(grid.scale);
grid.g=cell(grid.n_scale,1);
for i_scale = 1:parm.n_scale
    grid.g{i_scale}.sx = parm.scale(:,i_scale);
    grid.g{i_scale}.nx = 2.^grid.g{i_scale}.sx+1;
    grid.g{i_scale}.n = sum(grid.g{i_scale}.nx); % total number of cells
    grid.g{i_scale}.id=(1:grid.g{i_scale}.nx)';
    grid.g{i_scale}.dx = grid.grid.xmax ./ (grid.g{i_scale}.xn-1);

    for i_dim =1:grid.dim
        grid.g{i_scale}.x{i_dim}=linspace(0, grid.g.xmax(i_dim), grid.g{i_scale}.nx(i_dim))'; % coordinate of cells center
    end
    if grid.dim==2
        [grid.g{i_scale}.X{1}, grid.g{i_scale}.X{2}] = meshgrid(grid.g{i_scale}.x{1}, grid.g{i_scale}.x{2}); % matrix coordinate
    elseif grid.dim==3
        [grid.g{i_scale}.X{1}, grid.g{i_scale}.X{2}, grid.g{i_scale}.X{3}] = meshgrid(grid.g{i_scale}.x{1}, grid.g{i_scale}.x{2}, grid.g{i_scale}.x{3}); % matrix coordinate
    end
end
if ~isfield(parm, 'cstk_s') % cstk_s is the scale at which cst is switch on
    if ~isfield(parm, 'cstk'),      parm.cstk           = 1; end % constant path and kriging weight activated or not
    if parm.cstk
        parm.cstk_s = 0; % will never use cstk
    else
        parm.cstk_s = Inf; % will always use cstk
    end
end
if ~isfield(parm, 'fitvar'), parm.fitvar     =0; end % fit the variogram to the data or used the given one in parm.covar



% Kriging parameter
parm.k.covar = kriginginitiaite(parm.k.covar);
if ~isfield(parm, 'nscore'),        parm.nscore        =1; end % use normal score (strongly advice to use it.)
if ~isfield(parm, 'k') || ~isfield(parm.k, 'method'),  parm.k.method = 'sbss'; end
if ~isfield(parm, 'k') || ~isfield(parm.k, 'quad'),  parm.k.quad = 0; end
if ~isfield(parm, 'k') || ~isfield(parm.k, 'nb'),  parm.k.nb = [0 0 0 0 0; 5 5 5 5 5]; end
if ~parm.k.quad
    parm.k.n = sum(parm.k.nb(2,1:4));
end

if ~isfield(parm, 'k') || ~isfield(parm.k, 'sb') || ~isfield(parm.k.sb, 'nx') || ~isfield(parm.k.sb, 'ny') % super-block grid size (hard data)
    parm.k.sb.nx    = ceil(grid.g.x(end)/parm.k.covar(1).range(1)*3);
    parm.k.sb.ny    = ceil(grid.g.y(end)/parm.k.covar(1).range(2)*3);
end
if ~isfield(parm, 'k') || ~isfield(parm.k, 'wradius')
    parm.k.wradius  = Inf;
end
k = parm.k;

% Plot
if ~isfield(parm, 'plot') || ~isfield(parm.plot, 'bsgs'),  parm.plot.bsgs   =0; end
if ~isfield(parm, 'plot') || ~isfield(parm.plot, 'ns'),    parm.plot.ns     =0; end
if ~isfield(parm, 'plot') || ~isfield(parm.plot, 'sb'),    parm.plot.sb     =0; end
if ~isfield(parm, 'plot') || ~isfield(parm.plot, 'kernel'),parm.plot.kernel =0; end
if ~isfield(parm, 'plot') || ~isfield(parm.plot, 'fitvar'),parm.plot.fitvar =0; end
if ~isfield(parm, 'plot') || ~isfield(parm.plot, 'krig'),  parm.plot.krig   =0; end
if ~isfield(parm, 'plot') || ~isfield(parm.plot, 'presentation'),  parm.plot.presentation   =0; end



% Creation of the grid



%% * 1. *SUPERBLOCK GRID CREATION*
% A mask (Boolean value) of the hard data is assigned to each superblock
% as follow: Only the n-closest (normalized by the covariance range) points
% (inside the ellipse/windows) to the centre of the superblock will be
% true. During the kriging, the mask of the superblock of the estimated
% point will be used to select the hard to add to the kriging system
if strcmp(parm.k.method,'sbss') && Z.n > 0
    k.sb.nx = parm.k.sb.nx; % number of superblock grid
    k.sb.ny = parm.k.sb.ny;
    [k, Z] = SuperBlockGridCreation(k, grid.g.x(end), grid.g.y(end), Z, parm.plot.sb, parm.k.nb(2,5));
end



%% * 3. *NORMAL SCORE TRANSFORM*
% Based on the hard data (well samples), a normal score transform is
% created using interpolation with power or exponential tail extrapolation.
% Matlab symbolique function are used for efficient coding. The back
% transform of the prior normal distribution function is also created
% (return the pdf in the initial space from the mean and variance in the
% normal space)
Nscore = nscore({}, parm, parm.plot.ns);

% Create the normal space primary variable of known data
Z.d_ns = Nscore.forward(Z.d);



%% * 4. *RUN SIMULATION*

if strcmp(parm.k.method,'sbss'); k.sb.mask_ini = k.sb.mask; end % mask will be change after, we preserve its structure

if parm.par && parm.n_realisation~=1 % if parralelelisation is selected
    delete(gcp('nocreate'));
    poolobj=parpool(parm.par_n); % find the number of core available
    par_n_realisation = ceil(parm.n_realisation/poolobj.NumWorkers);
    
    RR=cell(poolobj.NumWorkers,1);
    tt=cell(poolobj.NumWorkers,1);
    
    parm_pool=parm;
    parm_pool.n_realisation = par_n_realisation;
    
    parfor pool_i=1:poolobj.NumWorkers
        [RR{pool_i}, tt{pool_i}]=SGS_in(Z, k, Nscore, grid.g, parm_pool);
    end
    delete(poolobj)
    
    Res=RR{1};
    for pool_i=2:numel(RR)
        for i_scale=1:parm.n_scale
            Res{i_scale}.m = [Res{i_scale}.m; RR{pool_i}{i_scale}.m];
            Res{i_scale}.m_ns = [Res{i_scale}.m_ns; RR{pool_i}{i_scale}.m_ns];
        end
        t.scale = [t.scale; tt{pool_i}.scale];
        t.cstk = [t.cstk; tt{pool_i}.cstk];
        t.pt = [t.pt; tt{pool_i}.pt];
        t.krig = [t.krig; tt{pool_i}.krig];
    end

else
    [Res,t_out] = SGS_in(Z, k, Nscore, grid.g, parm);
    t.scale = t_out.scale;
    t.cstk = t_out.cstk;
    t.pt = t_out.pt;
    t.krig = t_out.krig;
end


% save intial value
if strcmp(parm.k.method,'sbss'); k.sb.mask = k.sb.mask_ini; end



t.global = toc(t.global);



%% * 5. *SAVE IT*
filename=['result-SGSIM/', parm.familyname, 'SIM-', parm.name ,'_', datestr(now,'yyyy-mm-dd_HH-MM-SS'), '.mat'];
if parm.saveit
    mkdir(['result-SGSIM/', parm.familyname])
    save(filename, 'parm', 'Res', 'grid', 't', 'Prim', 'Prim_true', 'k',  'Nscore')
end

if parm.notify
    unix(['echo "Simulation ' parm.familyname ' has finish now (' datestr(now,'yyyy-mm-dd_HH-MM-SS') ') in ' num2str(t.global/60) 'min" | mail -s "Simulation Finish" ' parm.notify_email]);
    load handel; sound(y,Fs)
end

end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Res,t_out]=SGS_in(Z, k, Nscore, grid, parm)

clear parm.k

Res = cell(parm.n_scale,1);
t.toc.scale = [];
t.toc.cstk = [];
t.toc.pt = [];
t.toc.krig = [];

if parm.dist
    % k.covar(1).dist=[];
    Res{1}.dist=[];
    Res{1}.s =[];
end

if parm.plot.presentation
    i_frame=1;
    fig3=figure('position',[100 100 800 1000],'Color','white'); 
    %fig3=figure('position',[100 100 800 800],'Color','white'); 
end

if parm.varcovar % compute the empirical var-covariance matrix (Emery & Pelaez, 2011)
    % We store the lambda matrix in Res{end}, even for conditional simulation
    % we still store them all together on the same lambda (no lambda0).
    Res{end}.lambda = sparse(grid{parm.n_scale}.nxy,grid{parm.n_scale}.nxy);
    
    % But this require then to have all conditioning point in the finale grid!!
    assert(all(ismember(Z.x,grid{parm.n_scale}.x) & ismember(Z.y,grid{parm.n_scale}.y)), 'All conditioning data need to be on the last grid')
    
    % Assign to each hard data its index on the final grid
    Z.varcovar=zeros(Z.n,1);
    for i_hard_data=1:Z.n
        if ismember(Z.y(i_hard_data),grid{parm.n_scale}.y) && ismember(Z.x(i_hard_data),grid{parm.n_scale}.x) % if belong to the grid
            Z.varcovar(i_hard_data) = find(grid{parm.n_scale}.Y==Z.y(i_hard_data) & grid{parm.n_scale}.X==Z.x(i_hard_data));
        end
    end
    % We'll might need the index of hard data in the lambda matrix
    % (computing Covariance matrix.
    Res{end}.lambda_Prim = Z.varcovar;
end

%% * 1. *Simulation of Scale*
for i_scale=1:parm.n_scale % for each scale
    t.tic.scale = tic;
    
    %% * 1.1. *INITIATE SCALE SIMULATION*
    % Allocating space for resulting field.
    Res{i_scale}.x=grid{i_scale}.x;
    Res{i_scale}.y=grid{i_scale}.y;
    Res{i_scale}.X=grid{i_scale}.X;
    Res{i_scale}.Y=grid{i_scale}.Y;
    Res{i_scale}.nx=grid{i_scale}.nx;
    Res{i_scale}.ny=grid{i_scale}.ny;
    Res{i_scale}.nxy=grid{i_scale}.ny*Res{i_scale}.nx;
    Res{i_scale}.m=repmat({nan(grid{i_scale}.ny,grid{i_scale}.nx)},parm.n_realisation,1); % matrix des resutlats
    
    % Populate the grid from previous scale.
    if i_scale~=1 % not at first scale
        for i_sim=1:parm.n_realisation
            Res{i_scale}.m{i_sim}( 1:(grid{i_scale-1}.dy/grid{i_scale}.dy):end, 1:(grid{i_scale-1}.dx/grid{i_scale}.dx):end) = Res{i_scale-1}.m{i_sim};
        end
    end
    
    % Assimilate the hard data (Prim) into the grid
    hard_data_idx=find(ismember(Z.y,grid{i_scale}.Y)&ismember(Z.x,grid{i_scale}.X));
    for i_sim=1:parm.n_realisation
        for hd=1:numel(hard_data_idx)
            Res{i_scale}.m{i_sim}(Z.x(hard_data_idx(hd))==grid{i_scale}.X & Z.y(hard_data_idx(hd))==grid{i_scale}.Y) = Z.d(hard_data_idx(hd));
        end
    end
    
    % Should not remove from Prim. 
%     Prim.d(hard_data_idx)=[];
%     Prim.x(hard_data_idx)=[];
%     Prim.y(hard_data_idx)=[];
%     Prim.d_ns(hard_data_idx)=[];
    if parm.varcovar; Z.varcovar(hard_data_idx)=[]; end
    if strcmp(parm.k.method,'sbss'); k.sb.mask(:,:,hard_data_idx)=[]; end
    Z.n=numel(Z.d);
    
    % Create the normal space result matrix and populate with known value
    Res{i_scale}.m_ns = repmat({nan(size(Res{i_scale}.m{1}))},parm.n_realisation,1);
    I=~isnan(Res{i_scale}.m{1});
    for i_realisation=1:parm.n_realisation
        Res{i_scale}.m_ns{i_realisation}(I) = Nscore.forward(Res{i_scale}.m{i_realisation}(I));
    end
    
    % ... ?
    if parm.varcovar
        Res{i_scale}.varcovar_id=reshape( find(ismember(grid{parm.n_scale}.X, Res{i_scale}.x) &  ismember(grid{parm.n_scale}.Y, Res{i_scale}.y)), Res{i_scale}.ny, Res{i_scale}.nx);
    end
    
    
    
    %% * 1.2. *SPIRAL SEARCH*
    % Create the windows for kringing with the function to compute the
    % normalized distence and the order of visit of the cells. Spiral
    % Search setting: previously data (on grid{i_scale} location)
    if strcmp(parm.k.method,'sbss') %&& Prim.n>0
        k.ss.el.dw = ceil(min(max(k.covar(1).range*k.wradius./[grid{i_scale}.dx grid{i_scale}.dy]), max(grid{i_scale}.nx,grid{i_scale}.ny)));
        [k.ss.el.X, k.ss.el.Y] = meshgrid(-k.ss.el.dw:k.ss.el.dw);% grid{i_scale} of searching windows
        [k.ss.el.X_T, k.ss.el.Y_T]=rotredtrans(k.ss.el.X*grid{i_scale}.dx, k.ss.el.Y*grid{i_scale}.dy, k.covar(1).azimuth, k.covar(1).range); % transforms the grid{i_scale}
        k.ss.el.dist = sqrt(k.ss.el.X_T.^2 + k.ss.el.Y_T.^2); % find distence
        [k.ss.el.dist_s, k.ss.el.dist_idx] = sort(k.ss.el.dist(:)); % sort according distence.
        k.ss.el.dist_idx=k.ss.el.dist_idx(k.ss.el.dist_s<=k.wradius); % only keep those below the max radius;
        k.ss.el.dist_s=k.ss.el.dist_s(k.ss.el.dist_s<=k.wradius); % only keep those below the max radius;
        k.ss.el.X_s=k.ss.el.X(k.ss.el.dist_idx); % sort the axis
        k.ss.el.Y_s=k.ss.el.Y(k.ss.el.dist_idx);
        k.ss.el.n=numel(k.ss.el.X_s); %number of possible neigh
    end
    
    
    %% * 1.3. *Generate the order for visting cells*
    % Randomly permute the cell not known (to be visited). And generate the
    % Y.x and Y.y coordinate in a random order.
    
    if i_scale<parm.cstk_s
        parm.cstk=0;
    else
        parm.cstk=1;
    end
    
    [Res{i_scale}.sim] = definepath(Res{i_scale},grid{i_scale},parm);
    
    
    %% * 1.4. *RANDOM NORMAL FIELD*
    % This create the random Normal distribution used for sampling the
    % posteriori distribution at each point.
    % rng('shuffle')
    Res{end}.U = normcdf(randn(Res{i_scale}.sim.n,parm.n_realisation)); % random field
    
    
    %% * 2. *Simulation of Point*
    i_plot=0; % used for computing the number of point simulated. Used for ploting
    for i_pt=1:Res{i_scale}.sim.n; % loop over each point
        i_plot=i_plot+1;
        pt.i=i_pt;
        t.tic.pt = tic;
        
        if parm.cstk % if option constant weight is activate.
            % Find the current point position on the grid{i_scale}. It
            % will be used on each realisation.
            pt.x = Res{i_scale}.sim.x_r{1}(i_pt);
            pt.y = Res{i_scale}.sim.y_r{1}(i_pt);
            
            % Kriging system
            t.tic.krig = tic;
            pt.krig = kriging(pt,Res{i_scale},Z,k,parm,1);
            t.toc.krig = [t.toc.krig; toc(t.tic.krig)];
        end
        
        
        %% * 3. *Simulation of Realisation*
        for i_realisation=1:parm.n_realisation
            
            
            if ~parm.cstk
                % Find the current point position on the grid{i_scale}.
                % This is changing for each realisation
                pt.x = Res{i_scale}.sim.x_r{i_realisation}(i_pt);
                pt.y = Res{i_scale}.sim.y_r{i_realisation}(i_pt);
                
                % Kriging h_hist
                pt.krig = kriging(pt,Res{i_scale},Z,k,parm,i_realisation);
            end
            
            if ~isempty(pt.krig.mask.prim) || ~isempty(pt.krig.mask.res)
             pt.krig.m = pt.krig.lambda'* [Z.d_ns(pt.krig.mask.prim) ; Res{i_scale}.m_ns{i_realisation}(pt.krig.mask.res)];
            else % no neighbor
                pt.krig.m=0; % mean ?
            end
            
            % add the count to covar
            if parm.dist
                % sel_g_ini=[Prim.x Prim.y; Res{i_scale}.X(~isnan(Res{i_scale}.m{i_realisation}(:))) Res{i_scale}.Y(~isnan(Res{i_scale}.m{i_realisation}(:)))];
                % h_dist = pdist2([Res{i_scale}.x(pt.x) Res{i_scale}.y(pt.y)],sel_g_ini(pt.krig.mask,:))';
                % Res{1}.dist = [Res{1}.dist ; h_dist i_pt*ones(numel(h_dist),1)];
                Res{1}.s = [Res{1}.s;pt.krig.s];
            end
            
            %% * 3.1 *KRIGING DISTRIBUTION*
            % Back transform the normal distribution (from krigeage) in
            % original space in the  grid{i_scale}. this become
            % the prior distribution
            pt.krig.pdf = Nscore.dist(pt.krig.m, sqrt(pt.krig.s));
            pt.krig.pdf = pt.krig.pdf./sum(pt.krig.pdf);
            
            
            
            %% * 3.3  *SAMPLING*:
            % Sample a point in the posteri distribution. To do so the CDF
            % is created and interpolation tool is used to assure
            % randomized sampling.
            pt.krig.cdf = cumsum(pt.krig.pdf);
            
            % Interpolation only works if it is mono.. increasing. But the
            % cdf can reach 1 so we just add eps (very small value).
            if ~all(diff(pt.krig.cdf)>0)
                pt.krig.cdf = pt.krig.cdf +  linspace(0,numel(pt.krig.cdf)*eps*2,numel(pt.krig.cdf))';
                assert(all(diff(pt.krig.cdf)>0))
            end
            
            pt.sampled = interp1(pt.krig.cdf, Nscore.support_dist, Res{end}.U(i_pt,i_realisation),'pchip');
            
            
            %% * 3.4 *PLOTIT*
            if parm.plot.krig && i_realisation==1 && i_plot>1 %%  || Nscore.forward(pt.sampled)<-3  || Nscore.forward(pt.sampled)>3 ) %|| mod(i_plot,50)==0
                figure(1); clf
                
                subplot(3,2,[1 4]);
                hold on
                h1=imagesc(Res{i_scale}.x,Res{i_scale}.y,Res{i_scale}.m_ns{1},'AlphaData',~isnan(Res{i_scale}.m_ns{1}));
                
                if strcmp(parm.k.method,'sbss')
                    sb_i = min([round((Res{i_scale}.y(pt.y)-k.sb.y(1))/k.sb.dy +1)'; k.sb.ny]);
                    sb_j = min([round((Res{i_scale}.x(pt.x) -k.sb.x(1))/k.sb.dx +1)'; k.sb.nx]);
                    windows=nan(k.sb.ny,k.sb.nx);
                    for u=1:length(k.el_X_s) %.. look at all point...
                        for q=1:4 %... and assign it to the corresponding quadrant
                            if sb_i+k.qs(q,1)*k.el_X_s(u)<=k.sb.nx && sb_j+k.qs(q,2)*k.el_Y_s(u)<=k.sb.ny && sb_i+k.qs(q,1)*k.el_X_s(u)>=1 && sb_j+k.qs(q,2)*k.el_Y_s(u)>=1% check to be inside the grid
                                windows(sb_i+k.qs(q,1)*k.el_X_s(u), sb_j+k.qs(q,2)*k.el_Y_s(u))=1;
                            end
                        end
                    end
                    h2=imagesc(k.sb.x,k.sb.y,windows,'AlphaData',windows*.5);
                    h3=mesh([0 k.sb.x+k.sb.dx/2],[0 k.sb.y+k.sb.dy/2],zeros(k.sb.ny+1, k.sb.nx+1),'EdgeColor','k','facecolor','none');
                end
                
                tt=-pi:0.01:pi;
                x=Res{i_scale}.x(pt.x)+k.covar(1).range(1)*cos(tt);
                x2=Res{i_scale}.x(pt.x)+k.wradius*k.covar(1).range(1)*cos(tt);
                y=Res{i_scale}.y(pt.y)+k.covar(1).range(2)*sin(tt);
                y2=Res{i_scale}.y(pt.y)+k.wradius*k.covar(1).range(2)*sin(tt);
                h4=plot(x,y,'--r'); h5=plot(x2,y2,'-r');
                h6=plot([Res{i_scale}.x(pt.x) Res{i_scale}.x(pt.x)],[min(y2) max(y2)],'-r');
                h7=plot([min(x2) max(x2)], [Res{i_scale}.y(pt.y) Res{i_scale}.y(pt.y)],'-r');
                
                lambda_c= 36+60.*(abs(pt.krig.lambda)-min(abs(pt.krig.lambda)))./(range(abs(pt.krig.lambda))+eps); % eps avoid NaN
                
                h8=scatter(Z.x,Z.y,[],Z.d_ns,'s','filled');
                
                n_hd = numel(Z.x(pt.krig.mask.prim));
                sel=[Z.x(pt.krig.mask.prim) Z.y(pt.krig.mask.prim); Res{i_scale}.X(pt.krig.mask.res) Res{i_scale}.Y(pt.krig.mask.res)];
                XY_ns = [Z.d_ns(pt.krig.mask.prim) ; Res{i_scale}.m_ns{i_realisation}(pt.krig.mask.res)];
                h9=scatter(sel(1:n_hd,1),sel(1:n_hd,2),lambda_c(1:n_hd),XY_ns(1:n_hd),'s','filled','MarkerEdgeColor','k');
                h10=scatter(sel(n_hd+1:end,1),sel(n_hd+1:end,2),lambda_c(n_hd+1:end),XY_ns(n_hd+1:end),'o','filled','MarkerEdgeColor','k');
                h11=scatter(Res{i_scale}.x(pt.x),Res{i_scale}.y(pt.y),100,pt.krig.lambda'* XY_ns,'o','filled','MarkerEdgeColor','r','LineWidth',1.5);
                
                xlabel('x[m]');ylabel('y[m]');
                %colorbar;
                xlim([Res{i_scale}.x(1)-grid{i_scale}.dx/2 Res{i_scale}.x(end)+grid{i_scale}.dx/2])
                ylim([Res{i_scale}.y(1)-grid{i_scale}.dy/2 Res{i_scale}.y(end)+grid{i_scale}.dy/2])
                %set(gca,'YDir','reverse');
                
                % legend([h3 h6 h8 h9 h10 h11],'Super grid','Window search with quadrant','Hard data point','Selected hard data point','Selected previously simulated point','Simulated Point','Location','northoutside','Orientation','horizontal')
                
                
                subplot(3,2,5); hold on;
                plot( parm.support_dist,pt.krig.pdf)
                plot(  parm.support_dist, max(pt.krig.pdf)*pt.krig.cdf)
                plot([pt.sampled pt.sampled],[0 max(pt.krig.pdf)],'k')
                legend('Kiriging estimate','CDF','sampled')
                
                
                subplot(3,2,6); hold on;
                %[f,x]=hist(Prim_ini.d(:)); plot(x,f/trapz(x,f),'linewidth',2);
                %[f,x]=hist(Z.d(:)); plot(x,f/trapz(x,f),'linewidth',2);
                [f,x]=hist(Res{i_scale}.m{i_realisation}(:),50); plot(x,f/trapz(x,f));
                legend('Well sampling', 'ERT','Simulation')
                
                drawnow
                keyboard
            end
            
            if parm.plot.presentation
                c_axis=[-3 3];
                figure(fig3); clf;
                subplot(4,1,4); hold on;
                set(gca,'XTick',[],'YTick',[],'YColor',[38 50 56]/255,'XColor',[38 50 56]/255); box on;
                %xlabel('Z'); ylabel('PDF(Z)'); xlim(c_axis); ylim([0 .015]);caxis(c_axis); 
                xlim([0 Res{i_scale}.sim.n]);  ylim([0 1]); xlabel('Order of simulation'); ylabel('Time of Computation (sec.)')
                
                grid_width=0.5; % 1.5
                line_width = 2;
                line_color =[175 184 188]/255 ;%[207 216 220]/255
                node_size = 25; % 40
                node_full_size = 400; % 800
                
                subplot(4,1,[1 3]); 
                hold on;set(gca,'XTick',[]);set(gca,'YTick',[]);xlabel('X'); ylabel('Y'); axis off equal; caxis(c_axis)
                test = [255, 255, 204;161, 218, 180;65, 182, 196;44, 127, 184;8, 104, 172;37, 52, 148]/255;
                colormap([interp1(linspace(1,64,size(test,1)), test(:,1), 1:64)' interp1(linspace(1,64,size(test,1)), test(:,2), 1:64)' interp1(linspace(1,64,size(test,1)), test(:,3), 1:64)'])
                surf(grid{end}.X,grid{end}.Y,zeros(grid{end}.ny, grid{end}.nx),'EdgeColor',[207 216 220]/255,'facecolor','none','LineWidth',grid_width);
                plot(grid{end}.X(:),grid{end}.Y(:),'.','Color',[207 216 220]/255,'MarkerSize',node_size-5);axis equal tight
                scatter(grid{end}.X(~isnan(Res{i_scale}.m_ns{1})),grid{end}.Y(~isnan(Res{i_scale}.m_ns{1})),node_full_size,Res{i_scale}.m_ns{1}(~isnan(Res{i_scale}.m_ns{1})),'filled');
                plot(pt.x-1, pt.y-1,'o','MarkerSize',node_size,'Color',[0 145 234]/255,'LineWidth',2,'MarkerFaceColor','w');
                Res{end}.F(i_frame) = getframe(fig3);i_frame=i_frame+1;
                
                               
                % subplot(4,1,[1 3]);
                sel=[zeros(0,2) ; Z.x(pt.krig.mask.prim) Z.y(pt.krig.mask.prim); Res{i_scale}.X(pt.krig.mask.res) Res{i_scale}.Y(pt.krig.mask.res)];
                for i=1:size(sel,1)
                    plot([pt.x-1 sel(i,1)],[pt.y-1 sel(i,2)], '-','Color',line_color, 'linewidth',line_width);
                    plot(sel(i,1),sel(i,2), 'o', 'MarkerSize',node_size,'MarkerFaceColor','w', 'linewidth',2,'Color',line_color);
                end
                scatter(grid{end}.X(~isnan(Res{i_scale}.m_ns{1})),grid{end}.Y(~isnan(Res{i_scale}.m_ns{1})),node_full_size,Res{i_scale}.m_ns{1}(~isnan(Res{i_scale}.m_ns{1})),'filled');
                plot(pt.x-1, pt.y-1,'o','MarkerSize',node_size,'Color',[0 145 234]/255,'LineWidth',2,'MarkerFaceColor','w');
               
                % Res{end}.F(i_frame) = getframe(fig3);i_frame=i_frame+1;
                
                
                
                subplot(4,1,4);
                %plot( Nscore.support_dist, pt.krig.pdf,'Color',[0 145 234]/255,'LineWidth',1.5); % axis 'auto y'
                %Res{end}.F(i_frame) = getframe(fig3);i_frame=i_frame+1;

                %plot([pt.sampled pt.sampled],[0 interp1(Nscore.support_dist,pt.krig.pdf,pt.sampled)],'Color',[207 216 220]/255,'LineWidth',2)
                %scatter(pt.sampled,interp1(Nscore.support_dist,pt.krig.pdf,pt.sampled),[],pt.sampled,'o','filled')
                %Res{end}.F(i_frame) = getframe(fig3);i_frame=i_frame+1;
                
                hold on; plot(t.toc.pt,'-k'); plot(t.toc.pt,'.k'); 
                
                
                
                % subplot(4,1,[1 3]);
                % scatter(pt.x-1, pt.y-1,node_full_size,Nscore.forward(pt.sampled),'filled');
                Res{end}.F(i_frame) = getframe(fig3);i_frame=i_frame+1;

            end
            
            
            %% 3.5 Compute the empirical var-covariance matrix (Emery & Pelaez, 2011)
            if parm.varcovar && i_realisation==parm.n_realisation
                pt_id = find(Res{i_scale}.y(pt.y)==grid{parm.n_scale}.Y&Res{i_scale}.x(pt.x)==grid{parm.n_scale}.X);
                
                Res{end}.lambda(pt_id, [Z.varcovar(pt.krig.mask.prim); Res{i_scale}.varcovar_id(pt.krig.mask.res)]) = -pt.krig.lambda./sqrt(pt.krig.s);
                Res{end}.lambda(pt_id,pt_id) = 1/sqrt(pt.krig.s);
            end
            
            
            %% 3.6 *Back transform*
            Res{i_scale}.m{i_realisation}(pt.y,pt.x) = pt.sampled;
            Res{i_scale}.m_ns{i_realisation}(pt.y,pt.x) = Nscore.forward(pt.sampled); % add only the last simulated point to normal score
        end
        % NOTHING SHOULD BE DONE AFTER THAT WHICH INVOLVED THE WEIGHT
        % BECAUSE WE CHANGED RES.M...
        t.toc.pt=[t.toc.pt; toc(t.tic.pt)];
    end
    
    % display info
    % disp(['Simulation ' num2str(i_scale) '/' num2str(parm.n_scale) ' finished in : ' num2str(toc(t.tic.scale))])
    t.toc.scale = [t.toc.scale; toc(t.tic.scale)];
    
    
end
t_out = t.toc;

end
