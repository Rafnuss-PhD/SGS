function varargout = varcovar(Res,varargin)

assert( mod(nargin-1,2) == 0, 'need a pair number of input')
for i=1:2:numel(varargin)
    switch varargin{i}
        case 'CY'
            CY = varargin{i+1};
        case 'vario'
            vario = true(varargin{i+1});
            assert(islogical(vario),'vario is not a logical variable. 0 or 1 only')
        case 'ord'
            ord = true(varargin{i+1});
            assert(islogical(ord),'ord is not a logical variable. 0 or 1 only')
        case 'Prim'
            Prim = varargin{i+1};
        otherwise
            error('input not recognize')
    end
end

if ~exist('CY','var'); CY = 1; end;
if ~exist('vario','var'); vario = 0; end;
if ~exist('ord','var'); ord = 0; end;

% if presence of hard data, separate l from l0, compute EY.
if isfield(Res{end},'lambda_Prim') && numel(Res{end}.lambda_Prim)>0
    l = Res{end}.lambda;
    l(Res{end}.lambda_Prim,:) = [];
    l(:,Res{end}.lambda_Prim) = [];
    l0=full(Res{end}.lambda(:,Res{end}.lambda_Prim));
    l0(Res{end}.lambda_Prim,:) = [];
    
    if exist('Prim','var')
        EY = Res{end}.m{end};
        a=1:Res{end}.nxy;
        EY(~any(bsxfun(@eq,a,Res{end}.lambda_Prim),1)) = -inv(l) * l0 * Prim.d;
        varargout{4} = EY;
    end
    varargout{3} = l0;
else
    l=Res{end}.lambda;
end

% Compute CY
if CY==1
   CY = (l) \ transpose(inv(l));
elseif CY==0
   CY = [];
end



if vario
    % create the CY total matrix (CY_t) which also take the hard data into
    if isfield(Res{end},'lambda_Prim') % presence of hard data
        CY_t=zeros(size(CY)+numel(Res{end}.lambda_Prim));
        CY_t(Res{end}.lambda_Prim,:)=NaN; 
        CY_t(:,Res{end}.lambda_Prim)=NaN;
        CY_t(~isnan(CY_t)) = full(CY);
    else
        CY_t=CY;
    end

    
    % If ord, we remove the top triangular value of CY to get a ortiented
    % 2D variogram.
    if ord
        id=[];
        for i_scale=1:numel(Res)
            id = [id; Res{i_scale}.varcovar_id(Res{i_scale}.sim.xy_r)'];
        end
        
        % remove the top triangle of CY to have anisotropic variogram: CY_tri
        CY_ord = CY_t(id,id); clear CY_t
        [~,id2]=sort(id);
        CY_ord_tri = tril(CY_ord); clear CY_ord
        CY_ord_tri(CY_ord_tri==0) = NaN;
        CY_tri = CY_ord_tri(id2,id2); clear CY_ord_tri
        
        
        % Add back the hard data to get to total triangular CY
        if isfield(Res{end},'lambda_Prim') % presence of hard data
            CY_t=zeros(size(CY)+numel(Res{end}.lambda_Prim));
            CY_t(Res{end}.lambda_Prim,:)=NaN;
            CY_t(:,Res{end}.lambda_Prim)=NaN;
            CY_t(~isnan(CY_t)) = full(CY_tri);
        else
            CY_t=CY_tri;
        end
        clear CY_tri
        
    end
    
%     ny=Res{end}.ny;
%     nx=Res{end}.nx;
%     nxy=Res{end}.nxy;
% 
%     % create the 3D matrix of variogram. 
%     vario=cell(3,1);
%     vario{3} = nan(ny*2-1,nx*2-1,nxy);
%     for i=1:nxy
%         [id_y,id_x]=ind2sub([ny,nx],i);
%         vario{3}( (ny-id_y)+(1:ny) , (nx-id_x)+(1:nx),i) = reshape(full(CY_t(i,:)),ny,nx);
%     end
%     vario{2} = nanmean(vario{3},3);
% 
%  
%     vario{1}.h = [ vario{2}(ny,nx), nanmean( [vario{2}(ny,nx-1:-1:1); vario{2}(ny,nx+1:end)])];
%     vario{1}.v = [ vario{2}(ny,nx), nanmean( [vario{2}(ny-1:-1:1,nx)'; vario{2}(ny+1:end,nx)'])];

    D = pdist2([Res{end}.X(:) Res{end}.Y(:)],[Res{end}.X(:) Res{end}.Y(:)]);
    vario=struct;
    vario.h= unique(D);
    for i=1:numel(vario.h)
        id = D==vario.h(i);
        id2 = CY_t(id);
        vario.g_n(i) = sum(id(:));
        vario.g(i) = sum(id2)/vario.g_n(i);
        vario.std(i) = sqrt(  sum((id2(:)-vario.g(i)).^2) / (vario.g_n(i)-1)  );
    end
    
    varargout{5} = vario;
end


if ord
    id=[];
    for i_scale=1:numel(Res)
        id = [id; Res{i_scale}.varcovar_id(Res{i_scale}.sim.xy_r)'];
    end
    
    l = l(id,id);
    
    if CY
        CY = CY(id,id);
    end
    
    idx = triu(true(size(CY)), -1);
    CY(idx)=NaN;
    
end

varargout{1} = CY;
varargout{2} = l;




    