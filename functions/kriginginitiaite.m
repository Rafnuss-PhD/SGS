function covar = kriginginitiaite(covar)

% if isfield(covar, 'gen')% patch for old version of gen
%     parm=covar;
%     clear covar;
%     if isfield(parm.gen.covar, 'modele') && isfield(parm.gen.covar, 'c')
%         
%         covar=struct('model',num2cell(parm.gen.covar.modele(:,1)),...
%                 'c0',num2cell(parm.gen.covar.c)',...
%                 'range',num2cell(parm.gen.covar.modele(:,2:3)',1)'...
%                 ,'azimuth',num2cell(parm.gen.covar.modele(:,4)));
%     elseif isfield(parm.gen, 'covar')
%         % new version of gen
%         assert(isfield(parm.gen.covar, 'model'),'parm.covar is not properly define')
%         assert(isfield(parm.gen.covar, 'c0'),'parm.covar is not properly define')
%         assert(isfield(parm.gen.covar, 'range'),'parm.covar is not properly define')
%         assert(isfield(parm.gen.covar, 'azimuth'),'parm.covar is not properly define')
%         covar    = parm.gen.covar;
%         clear parm.gen
%     else
%         error('No parm.gen readable')
%     end
%     
% elseif isfield(parm, 'k') && isfield(parm.k, 'model') && isfield(parm.k, 'c')
%     % old version
%     covar=struct('model',num2cell(parm.k.model(:,1)),...
%         'c0',num2cell(parm.k.c)',...
%         'range',num2cell(parm.k.modele(:,2:3)',1)'...
%         ,'azimuth',num2cell(parm.k.model(:,4)));
% elseif isfield(parm.k, 'covar')
assert(isfield(covar, 'model'),'parm.covar is not properly define')
assert(isfield(covar, 'c0'),'parm.covar is not properly define')
assert(isfield(covar, 'range0'),'parm.covar is not properly define')
assert(isfield(covar, 'azimuth'),'parm.covar is not properly define')
% else
%     error('You need to define parm.k')
% end
% clear parm.k.model parm.k.c

for i=1:numel(covar)
    % Patch for old version
    switch covar(i).model
        case 1
            covar(i).model='nugget';
        case 2
            covar(i).model='exponential';
        case 3
            covar(i).model='gaussian';
        case 4
            covar(i).model='spherical';
    end

    switch covar(i).model
        case 'nugget'
            covar(i).g = @(h,r) h==0;
            intvario=1;
        case 'triangle'
            assert(numel(covar.range)==1,'only valid in 1D')
            intvario=1;
            covar(i).g = @(h) max(1-h,0);
        case 'circular'
            intvario=1.17;
            covar(i).g = @(h) 2/pi*(acos(min(h,1))-min(h,1).*sqrt(1-min(h,1).^2));
        case 'spherical'
             intvario=1.3;
            covar(i).g = @(h) 1-3/2*min(h,1)+1/2*min(h,1).^3;
        case 'cubic'
            intvario=1.43;
            covar(i).g = @(h) 1 - 7*min(h,1).^2 + 35/4*min(h,1).^3 - 7/2*min(h,1).^5 + 3/4*min(h,1).^7;
        case 'exponential'
            intvario = .41;
            covar(i).g = @(h) exp(-h);
        case 'gaussian'
            intvario=.58;
            covar(i).g = @(h) exp(-h.^2);
        case 'stable'
            assert(isfield(covar, 'alpha'),'alpha covar is not properly define')
            assert(covar(i).alpha>0 && covar(i).alpha<=2,'Alpha value not possible')
            intvario=.41;
            covar(i).g = @(h) exp(-(h).^covar(i).alpha);
        case 'power'
            assert(isfield(covar, 'alpha'),'alpha covar is not properly define')
            assert(covar(i).alpha>0 && covar(i).alpha<2,'Alpha value not possible')
            covar(i).g = @(h) 1-h.^covar(i).alpha;
            warning('Approx the integrale')
            intvario=1.5;
        case 'k-bessel'
            assert(isfield(covar, 'alpha'),'alpha covar is not properly define')
            assert(covar(i).alpha>0 && covar(i).alpha<2,'Alpha value not possible')
            intvario=[.35 .5];
            covar(i).g = @(h) 1/(2^(covar(i).alpha-1) * gamma(covar(i).alpha)) .* max(h,eps).^covar(i).alpha .* besselk(covar(i).alpha,max(h,eps));
        case 'logarithmic'
            covar(i).g = @(h) 1-log(h+1);
            warning('Approx the integrale')
            intvario=.7;
        case 'cauchy'
            assert(isfield(covar, 'alpha'),'alpha covar is not properly define')
            assert(covar(i).alpha>0,'Alpha value not possible')
            covar(i).g = @(h) (1+h.^2).^covar(i).alpha;
            warning('Approx the integrale')
            intvario=1;
        case 'hyperbolic'
            covar(i).g = @(h) 1./(1+h);
            warning('Approx the integrale')
            intvario=[.2 .05];
        case 'cardinal sine'
            intvario=.2;
            covar(i).g = @(h) sin(max(eps,h))./max(eps,h);
        case 'matheron'
            covar(i).g = @(h) 1./(h);
        otherwise
            error('Variogram type not defined')
    end
    
    if numel(covar(1).range0)==1 || numel(intvario)==1
        covar(i).range=covar(i).range0*intvario(1);
    elseif numel(covar(1).range0)==2
        covar(i).range=covar(i).range0*intvario(2);
    elseif numel(covar(i).range0)==3
        covar(i).range=covar(i).range0*intvario(3);
    end
    
end

for i=1:numel(covar)
    if numel(covar(1).range)==1 || numel(covar(1).azimuth)==0
        covar(i).cx = 1/diag(covar(i).range(1));
    elseif numel(covar(1).range)==2
        ang=covar(i).azimuth; cang=cos(ang/180*pi); sang=sin(ang/180*pi);
        rot = [cang,-sang;sang,cang];
        covar(i).cx = rot/diag(fliplr(covar(i).range));
    elseif numel(covar(i).azimuth)==3
        error('3D')
    end
end

end
