%% Accuracy of Chebyshev coefficients via aliasing
% Yuji Nakatsukasa, April 2016

%% 1. One dimension
% As shown in [1, Thm. 4.2], the Chebyshev coefficients
% of the degree $n$ polynomial interpolant
% $$ p=\sum_{i=0}^n \hat c_iT_i(x) $$
% of a (Lipschitz continuous) function $f$ with Chebyshev expansion 
% $$ f=\sum_{i=0}^\infty c_iT_i(x) $$
% are obtained by aliasing the coefficients for the 
% higher-degree Chebyshev polynomials. 
% Specifically, the coefficient $c_i$ has error 
% $$ \hat c_i- c_i = (c_{2n-i}+c_{4n-i}+\cdots)+(c_{2n+i}+c_{4n+i}+\cdots), $$
% except for the 0th coefficient,
% $$ \hat c_0-c_0 = c_{2n}+c_{4n}+\cdots, $$
% and the $n\mbox{th}$ coefficient, 
% $$ \hat c_n-c_n = c_{3n}+c_{5n}+\cdots. $$ 
 
%%
% Now, smooth functions have rapidly decaying
% Chebyshev coefficients, and in particular, if $f$ is analytic, 
% its Chebyshev coefficients decay geometrically at the rate
% $c_n=O(\mbox{exp}(-kn))$ for some constant $k>0$. 
% This geometric decay, combined with the aliasing formulae above,
% tells us that the Chebyshev coefficients of the interpolant
% $p$ have errors $\hat c_i-c_i$ of
% predictably varying magnitude, for the dominant term in the error is 
% $c_{2n}$ for $c_0-\hat c_0$, 
% $c_{3n}$ for $c_n-\hat c_n$, and $c_{2n-i}$ for the rest. 
% From this we expect that a Chebyshev interpolant
% for an analytic $f$ would have 
% high accuracy in $c_0$ and very high accuracy in $c_n$,
% relative to the other $c_i$. 

%% 
% Here is an illustration. 
% For an analytic function $f$, we find its Chebyshev expansion
% coefficients.
% We then compute a low-degree interpolant $p$
% and examine the accuracy of its Chebyshev coefficients. 

clear, close all
LW = 'linewidth'; MS = 'markersize'; FS = 'fontsize';
CO = 'color'; green = [0 .7 0];
lw = 2; ms = 10; fs = 16;

fori = @(x) log(sin(10*x)+2);

f = chebfun(fori);
p = chebfun(fori,round(length(f)/3));

pc = p.coeffs;
fc = f.coeffs;

plotcoeffs(f,'.',LW,lw,MS,ms), hold on
plotcoeffs(p,'r.',LW,lw,MS,ms)
plot(0:length(pc)-1,abs(pc-fc(1:length(pc)))+eps,'.',CO,green,LW,lw,MS,ms)

h_legend = legend('f','p','f-p');
set(h_legend,FS,fs)

%%
% The blue and red plots show the absolute values of the Chebyshev
% coefficients for $f$ and the interpolant $p$. 
% Our focus here is on the green dots, showing
% the error in Chebyshev coefficients.
% Two observations can be made:
% (i) the error grows geometrically with the degree until degree $n-1$, and 
% (ii) (note the rightmost green dot) at the end the error is 
% much much smaller.
% These effects clearly reflect the aliasing formulae given above. 

%{
% TO BE COMMENTED OUT OR DETAILED
% The high accuracy of $c_0$ is closely connected to Gauss quadrature. 
% If the
% function error $||f-p||$ is $\epsilon$, then $\hat c_0$ has 
% accuracy $O(\epsilon^2)$. Moreover, if $f$ is a polynomial of degree
% $2n-1$ or less, then $\hat c_0$ would be exactly equal to $c_0$,
% This is the same doubled-degree accuracy enjoyed by Gauss quadrature. 
%}

%% 
% The extremely high accuracy of the $n\mbox{th}$ coefficient 
% is peculiar to the Chebyshev interpolation process, and this
% coefficient is exact if 
% $f$ is a polynomial of degree up to $3n-1$.
% Although we are unaware of a practical application of this fact,
% it does suggest that if only the $n\mbox{th}$
% Chebyshev coefficient is of interest, then a $n+1$-point Chebyshev
% interpolation could suffice, even if the corresponding
% interpolant is rather poor. 

%%
% Let's repeat the computation with a non-analytic function.

fori = @(x)abs((x-0.1).^3); % twice differentiable but not analytic 

f = chebfun(fori);
p = chebfun(fori,round(length(f)/6));

pc = p.coeffs;
fc = f.coeffs;

clf
plotcoeffs(f,'.',LW,lw,MS,ms), hold on
plotcoeffs(p,'r.',LW,lw,MS,ms)
plot(0:length(pc)-1,abs(pc-fc(1:length(pc)))+eps,'.',CO,green,LW,lw,MS,ms)
xlim([0 length(f)/2])

h_legend = legend('f','p','f-p');
set(h_legend,FS,fs,'Location','Best')

%%
% For non-analytic functions the difference in accuracy is less
% prominent, because the Chebyshev coefficients decay more slowly. 
% The green plot is still roughly a 'mirrored' version of the blue
% one, and the first and $n\mbox{th}$ coefficients
% still have higher accuracy than most of the rest.
% Again, these are all consequences of the aliasing formula.

%% 2. Two dimensions 
% Let's try an analogous experiment in Chebfun2. As before we start with a
% smooth bivariate function, form its chebfun2, and construct
% a low-degree interpolant (degree 5 in each direction).
% Then we examine the accuracy in the coefficients,
% which we show in matrix form. 

p = chebfun2(@(x,y)sin(x+y)+cos(x-y));
pc = chebcoeffs2(p);

[x,y] = chebpts2(6); % construct chebfun2 of degree [5 5] from values on 6x6 grid 
pt = chebfun2(p(x,y));
ptc = chebcoeffs2(pt);

format shortE
ptc-pc(1:size(ptc,1),1:size(ptc,2))

%%
% The $(i,j)$ entry of this matrix is the error in the
% coefficient for $T_{j-1}(x)T_{i-1}(y)$. 
% Observations similar to the 1-d case can be made:
% looking horizontally or vertically, the low-degree coefficients are
% more accurate than the higher ones, except at the very ends, where the
% accuracy is much better. Note the exceptional accuracy at the lower-right
% corner. 

%% Reference
%
% [1] L. N. Trefethen, Approximation Theory and Approximation Practice,
% SIAM, 2013.
