%% Random functions in 2D
% Nick Trefethen, April 2015

%%
% (Chebfun example approx2/Random2D.m)
% [Tags: #randnfun2, #Chebfun2, #Chebfun3]

%%
% Recently Chebfun added the command |randnfun| for generating
% smooth random
% functions in 1D as well as the periodic |randfuntrig|.
% In keeping with Chebfun's mission of realizing continuous
% analogues of the discrete objects of Matlab, |randnfun| can
% be regarded as a continuous analogue of |randn|.
% Chebfun can construct 2D random functions too, with
% |randnfun2| and |randnfun2trig|.  Random functions
% in 3D or on the sphere or the disk have not yet been implemented.

%%
% A technical paper has not yete been written describing these
% functions, but in a word, the idea is that a ``smooth random function''
% is constructed from a finite Fourier series with independent
% normally distributed random coefficients. A parameter |dt| must
% be specified that specifies the associated space scale.  Approximately
% speaking, a random function contains wave numbers up to about
% $2\pi/dt$..

%%
% To illustrate, here is a random function with $dt = 0.2$ on
% a $2\times 1$ rectangle.  Negative values are black and
% positive values are white.
dt = 0.1;
rng(0), f = randnfun2(dt, [0 2 0 1]);
plot(f), view(0,90), colormap(gray(2))
axis equal, axis([0 2 0 1])

%%
% A contour plot shows more:
LW = 'linewidth';
contour(f,LW,2.5)
axis equal, axis([0 2 0 1])

%%
% To isolate the zero contours to high accuracy (though it takes
% longer), one could use |roots|.
c = roots(f);
plot(c,LW,3)
axis equal, axis([0 2 0 1])
axis equal, axis([0 2 0 1])

%%
% Here's a 3D plot
plot(f)

%% 
% Here for comparison is a periodic random function:
f = randnfun2trig(dt, [0 2 0 1]);
plot(f), view(0,90), colormap(gray(2))
axis equal, axis([0 2 0 1])

%%
% And here is a random function with $dt = 0.05$.
dt = 0.05;
f = randnfun2(dt, [0 2 0 1]);
plot(f), view(0,90), colormap(gray(2))
axis equal, axis([0 2 0 1])
