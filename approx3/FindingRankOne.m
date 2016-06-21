%% Finding a rank-one trivariate function basis
% Yuji Nakatsukasa, June 2016

%%
% (Chebfun example approx3/findingrankone.m)

%% 1. Rank-one trivariate functions 
% When a chebfun3 is constructed for a rank-one function 
% $f(x,y,z) = f_x(x)f_y(y)f_z(z)$, Chebfun is able to detect its numerical
% rank for efficient storage and subsequent computation.

f = chebfun3(@(x,y,z) sin(x).*cos(y).*exp(z));
rank(f)

%%
% A sum of $k(\geq 2)$ rank-one functions is usually of 
% rank (Tucker rank) $k$.  For example,
g = chebfun3(@(x,y,z) cos(x).*exp(y).*sin(z));
h = chebfun3(@(x,y,z) exp(x).*sin(y).*cos(z));
fhat = f+(g+h)/10;
rank(fhat)

%% 2. Finding a rank-one basis for subspace of trivariate functions
% Now consider the following problem.
% We are given the rank-3 function $\hat f$ above, 
% along with the rank-one functions $g$ and $h$, and 
% we would like to "recover" the rank-one function $f$
% (note that for rank-one functions, the Tucker and CP ranks are the same). 
% Note that $\hat
% f,g,h$ and $f,g,h$ span the same subspace of trivariate
% functions, so the goal can be rephrased as finding a rank-one basis for the 
% subspace spanned by $\hat f,g,h$. 
% This is a higher-order and continuous analogue of the problem considered in [1]. 
% A convenient way to obtain a rank-one function close to $\hat f$ is to do 
ftmp = chebfun3(@(x,y,z) fhat(x,y,z),'rank',[1 1 1]);

%%
% But this does not recover $f$ very well:
scale = f(1,1,1)/ftmp(1,1,1); % scalar scaling
plot(f-ftmp*scale),shg

%%
% Here is a simple algorithm, analogous to that in [2],
% that tries to recover the function $f$ more accurately. It is based on
% alternating projection between rank-one functions and 
% the subspace of trivariate functions spanned by $\hat f,g,h$.

n = length(f);
G = reshape(sample(g,n,n,n),[n^3,1]); % form vectors of values at Chebyshev tensor grid
H = reshape(sample(h,n,n,n),[n^3,1]);
F = reshape(sample(fhat,n,n,n),[n^3,1]);
[Q,~] = qr([G H F],0);                % Q is the subspace spanned by $fhat,g,h$

for ii = 1:10
    Ftmp = reshape(sample(ftmp,n,n,n),[n^3,1]);
    Ftmp = Q*(Q'*Ftmp);               % projection onto subspace
    ftmp = chebfun3(reshape(Ftmp,n,n,n),'rank',[1 1 1]); % proj onto rank-1 funs     
    scale = f(1,1,1)/ftmp(1,1,1);     % scalar scaling
    err(ii) = norm(f-ftmp*scale);
end

plot(f-ftmp*scale),shg

%%
% The last plot suggests convergence of ftmp to $f$. 
% Indeed, under mild assumptions and with an initial
% guess close to an intersection
% point, alternating directions converges linearly to the intersection. 
% For this example; the
% convergence of $\|f-\hat f\|$ is convincingly linear.
MS = 'Markersize'; ms = 18;LW = 'linewidth'; 
LW = 'linewidth'; MS = 'markersize'; FS = 'fontsize';
TEX = 'interpreter';tex = 'latex'; 
lw = 2; ms = 12; fs = 14; ffs = 12;

semilogy(err,'-o',LW,lw)
xlabel('iteration',FS,fs)
ylabel('error $\|f-\hat f\|$',TEX,tex,FS,fs)

%% 3. References
%
% [1] D. Drusvyatskiy, A. D. Ioffe, and A. D. Lewis, 
% Transversality and alternating projections
% for nonconvex sets, _Found. Comput. Math._ 15 (2015), 1637-1651.
% 
% [2] Y. Nakatsukasa, T. Soma, and A. Uschmajew, Finding a low-rank
% basis in a matrix subspace, _Mathematical Programming_, to appear.
