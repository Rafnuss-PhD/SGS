function [Res, t, k, parm, filename] = SGS_par(Z, parm)

addpath('./functions/.')
if ismac
    % Code to run on Mac plaform
elseif isunix
    addpath('./functions/MinMaxSelection-UNIX/.')
elseif ispc
    addpath('./functions/MinMaxSelection-PC/.')
    addpath('./functions/partialSort-PC/.')
else
    disp('Platform not supported')
end
t.global = tic;
% addpath(genpath('./.'))

%% * *INPUT CEHCKING*
% This section of the code generates a valid parm structure based on the
% inputted parm. If some parm value haven't been input, this section will
% fill automatically with defautl value. This may not allowed be the best.

% Force input in column
Z.x=Z.x(:);Z.y=Z.y(:);Z.d=Z.d(:);

% Paramter settings
if ~isfield(parm, 'seed'),          parm.seed           = 'shuffle'; end
rng(parm.seed);
if ~isfield(parm, 'saveit'),        parm.saveit         = 1; end % bolean, save or not the result of simulation
if ~isfield(parm, 'name'),          parm.name           = ''; end % name use for saving file
if ~isfield(parm, 'familyname'),    parm.familyname     = ''; end
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
if ~isfield(parm, 'varcovar'),      parm.varcovar        = 0; end

% Scale and weight parameters
if ~isfield(parm, 'nx') || ~isfield(parm.g, 'scale')
    parm.g.scale=[6 ; 6];
end
if ~isfield(parm, 'g') || ~isfield(parm.g, 'max')
    parm.g.max=[2^parm.g.scale(1,end) 2^parm.g.scale(2,end)];
end

if ~isfield(parm, 'cstk_s') % cstk_s is the scale at which cst is switch on
    if ~isfield(parm, 'cstk'),      parm.cstk           = 1; end % constant path and kriging weight activated or not
    if parm.cstk
        parm.cstk_s = 0; % will never use cstk
    else
        parm.cstk_s = Inf; % will always use cstk
    end
end

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
    parm.k.sb.nx    = ceil( parm.g.max(1)/parm.k.covar(1).range(1)*3);
    parm.k.sb.ny    = ceil( parm.g.max(2)/parm.k.covar(1).range(2)*3);
end
if ~isfield(parm, 'k') || ~isfield(parm.k, 'wradius')
    parm.k.wradius  = Inf;
end
k = parm.k;


% Check the input for correct size dans dimension
assert(size(Z.d,2)<=1,'X.d is not a vertical vector 1D');
assert(size(Z.x,2)<=1,'X.dx is not a vertical vector 1D');
assert(size(Z.y,2)<=1,'X.y is not a vertical vector 1D');
assert(all(size(Z.y)==size(Z.x)),'X.x and X.y don''t have the same dimension');
assert(all(size(Z.d)==size(Z.x)),'X.d and X.x (or X.y) don''t have the same dimension');
assert(isfield(parm, 'nx'), 'Specify a number of cell in x' )
assert(isfield(parm, 'ny'), 'Specify a number of cell in y' )

% Creation of the grid an path
parm.x=1:nx;
parm.y=1:ny;
[parm.X, parm.Y] = meshgrid(parm.x,parm.y);

[Res{i_scale}.sim] = definepath(parm);




%% * 1. *SUPERBLOCK GRID CREATION*
% A mask (Boolean value) of the hard data is assigned to each superblock
% as follow: Only the n-closest (normalized by the covariance range) points
% (inside the ellipse/windows) to the centre of the superblock will be
% true. During the kriging, the mask of the superblock of the estimated
% point will be used to select the hard to add to the kriging system
if strcmp(parm.k.method,'sbss') && Z.n>0
    k.sb.nx = parm.k.sb.nx; % number of superblock grid
    k.sb.ny = parm.k.sb.ny;
    [k, Z] = SuperBlockGridCreation(k, parm.g.max(1), parm.g.max(2), Z, parm.plot.sb, parm.k.nb(2,5));
end


end
