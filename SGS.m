%% Sequential Gaussian Simulation
% Traditional SGS produced Gaussian realization 
%
% *INPUT* 
%
% * nx,ny: number of cell in x and y [int]
% * m: number of realization [int]
% * covar: covariance function structure
% * neigh: neighborhood search strategie
% * parm: Other paramerter

% 
% *OUPUT*
%
% * Rest: Realizations matrix (nx,ny,m)
% * t: time structure
%%

function [Rest, t] = SGS(nx,ny,m,covar,neigh,parm, method,hd)

addpath('./functions/')
 
%% 1. Checking and formating Input
% This section of the code generates a valid parm structure based on the
% inputted parm. If some parm value haven't been input, this section will
% fill automatically with defautl value. This may not allowed be the best.

% Grid size and realization number
assert(floor(nx)==nx && nx>0,'nx must be a positive integer');
assert(floor(ny)==ny && ny>0,'nx must be a positive integer');
assert(floor(m)==m && m>0,'The number of realizations must be a positive integer');

% Covariance structure. See |kriginginitiaite.m| for more information
if ~isfield(covar, 'model') ,  covar.model = 'spherical'; end
if ~isfield(covar, 'c0') ,  covar.c0 = 1; end
if ~isfield(covar, 'range0') ,  covar.range0 = [nx/5 ny/5]; end
if ~isfield(covar, 'azimuth') ,  covar.azimuth = 0; end
if ~isfield(covar, 'alpha') ,  covar.alpha = 1; end
covar = kriginginitiaite(covar);

% Neighborhood search strategie
if ~isfield(neigh, 'method'),  neigh.method = 'sbss'; end
if ~isfield(neigh, 'lookup'),  neigh.lookup = false; end
if ~isfield(neigh, 'nb'),  neigh.nb = 30; end
if ~isfield(neigh, 'wradius') neigh.wradius  = 3; end

% Paramter settings
if ~isfield(parm, 'seed_path'),     parm.seed_path      = 'shuffle'; end
if ~isfield(parm, 'seed_U'),        parm.seed_U         = 'shuffle'; end
if ~isfield(parm, 'seed_search'),   parm.seed_search         = 'shuffle'; end
if ~isfield(parm, 'saveit'),        parm.saveit         = 0; end % bolean, save or not the result of simulation
if ~isfield(parm, 'name'),          parm.name           = ''; end % name use for saving file
if ~isfield(parm, 'path'),          parm.path            = 'linear'; end
if ~isfield(parm, 'path_random'),   parm.path_random     = 1; end
if ~isfield(parm, 'mg'),            parm.mg              = 1; end


%% 2. Use the requested method or find the optimal one

switch method
    case 'trad'
        [Rest, t] = SGS_trad(nx,ny,m,covar,neigh,parm);
    case 'cst'
        [Rest, t] = SGS_cst(nx,ny,m,covar,neigh,parm);
    case 'cst_par'
        [Rest, t] = SGS_cst_par(nx,ny,m,covar,neigh,parm);
    case 'cst_par_cond'
        hd.x=hd.x(:); hd.y=hd.y(:); hd.d=hd.d(:);
        if ~isfield(hd, 'n'), hd.n= numel(hd.d); end
        if ~isfield(hd, 'id')
            hd.id = sub2ind([ny nx],hd.y,hd.x);
        else
            hd.id=hd.id(:);
        end
        hd.scale=nan(hd.n,1);
        [Rest, t] = SGS_cst_par_cond(nx,ny,m,covar,neigh,parm,hd);
    case 'hybrid'
        if ~isfield(parm, 'hybrid'), parm.hybrid=5;end
        [Rest, t] = SGS_hybrid(nx,ny,m,covar,neigh,parm);
    case 'varcovar'
        Rest = SGS_varcovar(nx,ny,m,covar,neigh,parm);
    otherwise
        
end


%% 3. Saving if requested
if parm.saveit
    filename=['result-SGS/SIM-', parm.name ,'_', datestr(now,'yyyy-mm-dd_HH-MM-SS'), '.mat'];
    mkdir('result-SGS/')
    save(filename, 'parm','nx','ny', 'Rest', 't', 'k','U')
end
