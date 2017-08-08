function [pt_krig] = kriging(pt,Res,Prim,krig,parm,i_realisation)
%% kriging return the mean and variance estimate at the point Y.pt
%
% INPUT:
% * Y       : Primary variable
% * X       : Secondary variable
% * k       : kriging information
%
% OUTPUT:
% * krig_m  : kriging mean estimate
% * krig_s  : kriging variance estimate
%
% * *Author:* Raphael Nussbaumer (raphael.nussbaumer@unil.ch)
% * *Date:* 02.02.2015

% *SELECTION OF NEIGHBOURING POINT*
pt_krig = neighsearch(pt,Res,Prim,krig,parm,i_realisation);

sel=[zeros(0,2) ; Prim.x(pt_krig.mask.prim) Prim.y(pt_krig.mask.prim); Res.X(pt_krig.mask.res) Res.Y(pt_krig.mask.res)];
    
% * *KRIGING*: Find his kringing value in noraml space:

a0_C=covardm_perso(sel,[Res.x(pt.x) Res.y(pt.y)],krig.covar);
ab_C=covardm_perso(sel,sel,krig.covar);

pt_krig.lambda = ab_C \ a0_C;
pt_krig.s = sum([krig.covar.c0]) - pt_krig.lambda'*a0_C;

% disable for time saving
 assert(~any(isnan(pt_krig.lambda)),'the kriging coeff is NaN')
 assert(pt_krig.s>0,'the kriging std result is less than zero')

end

    