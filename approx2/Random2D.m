%% Random functions in 2D
% Nick Trefethen, April 2017

%%
% (Chebfun example approx2/Random2D.m)
% [Tags: #randnfun2, #Chebfun2, #Chebfun3]

%%
% Recently Chebfun added the command |randnfun| for generating
% smooth random
% functions in 1D, as well as its periodic cousin |randfuntrig|.
% In keeping with Chebfun's mission of realizing continuous
% analogues of the familiar discrete objects, |randnfun| can
% be regarded as a continuous analogue of the Matlab command |randn|.
% Chebfun can construct 2D random functions too, with
% |randnfun2| and |randnfun2trig|.  Random functions
% in 3D or on the sphere or the disk have not yet been implemented.

%%
% A technical paper has not yet been written describing these
% functions, but in a word, the idea is that a "smooth random function"
% is constructed from a finite Fourier series with independent
% normally distributed random coefficients. A parameter dt must
% be specified that sets the associated space scale.  Approximately
% speaking, a random function contains wave numbers up to about
% $2\pi/dt$.

%%
% To illustrate, here is a random function with $dt = 0.4$ on
% a $2\times 1$ rectangle.  Negative values are black and
% positive values are white.
dt = 0.2;
rng(0), f = randnfun2(dt, [0 2 0 1]);
plot(f), view(0,90), colormap(gray(2))
axis equal, axis([0 2 0 1])
LW = 'linewidth'; XT = 'xtick'; YT = 'ytick';
set(gca,XT,0:.5:2,YT,0:.5:1)

%%
% A contour plot shows more:
contour(f,LW,2.5), colormap('default'), colorbar
axis equal, axis([0 2 0 1])
set(gca,XT,0:.5:2,YT,0:.5:1)

%%
% To isolate the zero contours to high accuracy (though it takes
% longer), one could use |roots|.
c = roots(f);
plot(c,LW,3)
axis equal, axis([0 2 0 1])
set(gca,XT,0:.5:2,YT,0:.5:1)

%%
% Here's a 3D plot.
plot(f)
view(-20,50), camlight left

%% 
% Here for comparison is a periodic random function.
f = randnfun2trig(dt, [0 2 0 1]);
plot(f), view(0,90), colormap(gray(2))
axis equal, axis([0 2 0 1])
set(gca,XT,0:.5:2,YT,0:.5:1)

%%
% And here are random functions with $dt = 0.1$
dt = 0.1; f = randnfun2(dt, [0 2 0 1]);
plot(f), view(0,90), colormap(gray(2))
axis equal, axis([0 2 0 1])
set(gca,XT,0:.5:2,YT,0:.5:1)

%%
% and with $dt = 0.05$
dt = 0.05; f = randnfun2(dt, [0 2 0 1]);
plot(f), view(0,90), colormap(gray(2))
axis equal, axis([0 2 0 1])
set(gca,XT,0:.5:2,YT,0:.5:1)
