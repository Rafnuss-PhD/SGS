%% NSCORE_PERSO
% This function is computing the Normal zscore transform of the input
% vector the function return two fonction handle : one for the normal transform
% and the other one for the back-transform. The function use the inputed
% vector to create the normal transform. Then using interpolation, the
% function handle are created.
%
% INPUT:
% * X            : Input vector
% * support_dist : Vector of value on which the distribution is built
%
% OUTPUT:
% * NscoreT      :
% * NscoreTinv   :
%
% * *Author:* Raphael Nussbaumer (raphael.nussbaumer@unil.ch)
% * *Date:* 29.01.2015

function Nscore = nscore(kern, parm, plotit) % X,support_dist,method,extrapolationMethod,plotit



if parm.nscore
    prior = kern.prior./sum(kern.prior);
    prior(prior<eps)=2*eps;
    cdf = cumsum(prior) ./sum(prior);
    
    Nscore.T_F = griddedInterpolant(kern.axis_prim,cdf,'pchip','pchip');
    Nscore.Tinv_F = griddedInterpolant(cdf,kern.axis_prim,'pchip','pchip');

    Nscore.inverse = @(y) Nscore.Tinv_F(normcdf(y)); % back-transform a value in normal space by taking the normcdf.
    Nscore.forward = @(y) norminv( Nscore.T_F(y) );
    
    kernel_y_ns = norminv(cdf);
    kernel_y_ns_mid = ( kernel_y_ns(1:end-1)+kernel_y_ns(2:end) ) /2;
    Nscore.dist = @(mu,sigma) [ normcdf(kernel_y_ns_mid,mu,sigma) ; 1] - [0 ; normcdf(kernel_y_ns_mid(),mu,sigma)]; 

else
    Nscore.support_dist = linspace(-5,5,500)';
    Nscore.forward = @(x) x';
    Nscore.inverse = @(x) x';
    Nscore.dist    = @(mu,sigma) normpdf(Nscore.support_dist,mu,sigma)/sum(normpdf(Nscore.support_dist,mu,sigma));
end





return
%% NOTE:
% Here are 6 Method to compute the transform and back transform
% A:  input vector
% B: Normal z-score of A
% b: point of a std normal distribution
% a: the back transform of b


% Method 1: Inital script from Paolo
B=nscoretool(A);
nt=length(A);
zmin=min(b);
zmax=max(b);
ltail=2;
ltpar=2;
utail=1;
utpar=2;
a=backtrtool_pr(b,nt,A,B,zmin,zmax,ltail,ltpar,utail,utpar);

plot(A,B,'o')



% Method 2: mGstat toolbox
% w1,dmin : Extrapolation options for lower tail. w1=1 -> linear interpolation, w1>1 -> Gradual power interpolation
% w2,dmax : Extrapolation options for lower tail. w1=1 -> linear interpolation, w1<1 -> Gradual power interpolation
% DoPlot : plot
% d_nscore : normal score transform of input data
% o_nscore : normal socre object containing information needed to perform normal score backtransform.
[B,o_nscore]=nscore(A,w1,w2,dmin,dmax);
a=inscore(b,o_nscore);



% Method 3. ECDF
CDF_inv=norminv(ecdf(A));
B = CDF_inv(tiedrank(A));
a=quantile(A,normcdf(b));



% Method 4. TieRank
B = norminv( tiedrank(A)/(numel(A)+1));
a=quantile(A,normcdf(b));



% Method 5. http://ch.mathworks.com/help/stats/examples/nonparametric-estimates-of-cumulative-distribution-functions-and-their-inverses.html#zmw57dd0e1074
[Bi,xi] = ecdf(A);
n=numel(A);
xj = xi(2:end);
Bj = (Bi(1:end-1)+Bi(2:end))/2;
xj = [xj(1)-Bj(1)*(xj(2)-xj(1))/((Bj(2)-Bj(1)));  xj;  xj(n)+(1-Bj(n))*((xj(n)-xj(n-1))/(Bj(n)-Bj(n-1)))];
Bj = [0; Bj; 1];

F    = @(y) norminv(interp1(xj,Bj,y,'linear','extrap'));
Finv = @(u) normcdf(interp1(Bj,xj,u,'linear','extrap'));

B=norminv(F(A));
a=Finv(normcdf(b));




% Method 6. http://ch.mathworks.com/help/stats/examples/nonparametric-estimates-of-cumulative-distribution-functions-and-their-inverses.html#zmw57dd0e1147
Fy = ksdensity(A, A, 'function','cdf', 'width',.35);
Finv = @(u) interp1(Fy,y,u,'linear','extrap');

B=norminv(Fy);
a=Finv(normcdf(b));


% Method 7 Complexe
% Compute the empirical cdf
[ft,xt] = ecdf(X);

% ecdf compute the step function because it is empirical, and there is
% one step per value of X (the height of the step is constent). The real
% pdf of this variable is something (no exacly, but closer) to a linear
% extrapolation between the center of the step. In order to get these
% point, we average the value of x and f:
x = xt(2:end);
f = (ft(1:end-1)+ft(2:end))/2;
% figure; ecdf(X); hold on; plot(x,f,'o-r')


% Now, our problem is to find a good extrapolation at each tail. This is
% important because of that : the prior pdf is estimate in the normal space and back
% transform to initial space to do the sampling. Then the value is
% transform to to normal space for the next prior estimate (kriging).
%
% After trying this:
% x = [support_dist(1)-1 ; x(1)-f(1)*(x(2)-x(1))/((f(2)-f(1)));  x;  x(end)+(1-f(end))*((x(end)-x(end-1))/(f(end)-f(end-1))) ; support_dist(end)+1];
% f = [0+eps; 0+2*eps; f; 1-2*eps; 1-eps];
% x = [support_dist(1); x; support_dist(end)];
%f = [eps; f; 1-eps];
%
% The solution was to fit a expenential, between the first point of the
% kernel and the first point of the ecdf, then use the point of the kernel
% as ...




dx=min([10; numel(X)]); % number of point to take in the interpolation:
if support_dist(1)>0
    fx= 'power1';
elseif support_dist(1)==0
    support_dist(1)=eps;
    fx= 'power1';
else
    fx='exp1';
end

% Lower Tail
w=ones(dx+1,1); w(2)=1000000; % set a very high weight on the last poin to force to be exaclty there (needed for interpolation later)
options = fitoptions('Weights',w); % trick to force the point at x(1)
f_low=fit([support_dist(1); x(1:dx)],[eps; f(1:dx)], fx,options);
% figure;plot(f_low,[support_dist(1); x(1:dx)],[eps; f(1:dx)])

% Upper Tail
w=ones(dx+1,1); w(end-1)=1000000;
options = fitoptions('Weights',w);
f_high=fit([x(end-dx+1:end); support_dist(end)],1-[f(end-dx+1:end); 1-eps],fx,options);
% figure; plot(f_high,[x(end-dx-1:end); support_dist(end)],1-[f(end-dx-1:end); 1-eps])

% Add point upper and lower tail for the interpolation. The point added are
% at support_dist.
x2 = [support_dist(support_dist<x(1))   ; x ; support_dist(support_dist>x(end))];

f2_b=feval(f_low,support_dist(support_dist<x(1)));
f2_e=1-feval(f_high,support_dist(support_dist>x(end)));


f2 = [ f2_b  ; f ; f2_e ];

% check for monoticity in order to interpolate later
i=0;
while ~all(diff(f2)>0)
    id=find(diff(f2)<0);
    f2(id+1) = f2(id+1) + eps;
    i=i+1;
    if i>10
        error('too many correction...')
    end
    f2(id+[-2:2])
end

% f2 and x2 are the vector on which the interpolation between data and
% their probabilite is made.
% figure; ecdf(X); hold on; plot(x,f,'o'); plot(x2,f2,'-d')

Nscore.T_F = griddedInterpolant(x2,f2,method,extrapolationMethod);
Nscore.Tinv_F = griddedInterpolant(f2,x2,method,extrapolationMethod);

% Function of Normal Score Transform. 
% We need to associate to each probability (from the ecdf) a value in the
% normal space which correspond.
Nscore.inverse = @(y) Nscore.Tinv_F(normcdf(y)); % back-transform a value in normal space by taking the normcdf.
% Nscore.forward = @(y) norminv( min([ 1-eps*ones(numel(y),1) max([ eps*ones(numel(y),1) Nscore.T_F(y)],[],2) ],[],2));
Nscore.forward = @(y) norminv( Nscore.T_F(y) );


% Kriging generate a mean and standard deviation in the std normal space.
% We want to transform this std normal distribution in the original space.
% The final distribution is on the grid of support_dist. Therefore we compute
% the nscore transform of the support_dist cell and compute the probability
% corresponding of the normal distribution generated by kriging (mu, sigma)
kernel_y_ns = Nscore.forward(support_dist);

% Method with normpdf
% kernel_b = (kernel_y_ns(3:end)-kernel_y_ns(1:end-2) )/2;
% kernel_b = [kernel_b(1);kernel_b;kernel_b(end)];
% Nscore.dist = @(mu,sigma) normpdf(kernel_y_ns,mu,sigma).*kernel_b;

% Method with normcdf
kernel_y_ns_mid = ( kernel_y_ns(1:end-1)+kernel_y_ns(2:end) ) /2;
Nscore.dist = @(mu,sigma) [ normcdf(kernel_y_ns_mid,mu,sigma) ; 1] - [0 ; normcdf(kernel_y_ns_mid(),mu,sigma)];


end
