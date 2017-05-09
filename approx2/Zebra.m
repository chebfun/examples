%% Zebra plots
% Nick Trefethen, May 2017

%%
% (Chebfun example approx2/Zebra.m)

%%
% Instead of a plot showing
% many function values, we may prefer to show just a
% plus/minus distinction.  For this there is the
% |'zebra'| option in Chebfun2, Spherefun, and Diskfun.

%%
% For example, here is zebra plot of a certain function on
% the disk.
% For fun we've changed the colors from the usual black/white.
cheb.xydisk;
f = sin(20*(x+y).*(1+y));
plot(f, 'zebra')
colormap([1 1 0; 0 0 0])
axis off

%%
% Normally, however, the plots show zebras rather than bumblebees.
% Negative values are black and positive values are white.
% Here is an example on the sphere.
f = spherefun.sphharm(15,5);
plot(f,'zebra')
axis off

%%
% And here is an example on the unit square.
f = randnfun2(.1);
plot(f, 'zebra')
axis equal off
