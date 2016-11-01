%% Partitioning spherefuns
% Behnam Hashemi, November 2016

%%
% (Chebfun Example sphere/SpherefunPartition.m)
% [Tags: #Spherefun]

%%
% Assume that $f(x,y,z)$ is a function defined over the unit 2-sphere
% in three dimensions. Our aim is to explore the building blocks of $f$ 
% using the _partition_ command. Let's start with a spherefun object:
f = spherefun(@(x,y,z) sinh(5*x.*y.*z).*cos(x-y+2*z))
plot(f), axis off, hold on
contour(f, 'color','k'), 

%%
% A spherefun can be seen as a sum of two spherefuns, one of them 
% even/$\pi$-periodic and the other odd/$\pi$-anti-periodic [1]. Recall that
% a univariate function $g$ is $\pi$-anti-periodic if $g(x+\pi) = -g(x)$. 
% The command `[fep, foa] = partition(f)' partitions $f$ accordingly.
[fep, foa] = partition(f)
err = norm(fep+foa - f)

%%
% _fep_ has a CDR decomposition [1] whose columns are even and whose rows 
% are $\pi$ periodic (not just $2\pi$!):
[Ce, D, Rp] = cdr(fep);
clf, plot(Ce)
grid on, title('Columns of the even part of f')

%%
clf, plot(Rp)
grid on, title('Rows of the \pi-periodic part of f')

%%
% The other part of $f$, _foa_, has a CDR decomposition whose columns are 
% odd and whose rows are $\pi$ anti-periodic:
[Co, D, Ra] = cdr(foa);
plot(Co), 
grid on, title('Columns of the odd part of f')

%%
clf, plot(Ra)
grid on, title('Rows of the \pi-antiperiodic part of f')

%% References
% 1. A. Townsend, H. Wilber, and G. Wright, Computing with functions in 
% spherical and polar geometries I. The sphere. _SIAM J. Sci. Comput._,
% 38 (2016) C403-C425.