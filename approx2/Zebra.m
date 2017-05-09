%% Zebra plots
% Nick Trefethen, May 2017

%%
% (Chebfun example approx2/Zebra.m)

%%
% Sometimes less is more.  Instead of a plot showing
% many function values, we may prefer to show just a
% plus/minus distinction.  For this there is the
% |'zebra'| option in Chebfun2, Spherefun, and Diskfun.

%%
% For example, here is zebra plot of a certain function on
% the disk.  Negative values are black and positive values are white.
cheb.xydisk;
f = sin(20*(x+y).*(1+y));
plot(f, 'zebra')
axis off

%%
% Here is an example on the sphere.
f = spherefun.sphharm(15,5);
plot(f,'zebra')
axis off

%%
% And here is a rectangle.
f = randnfun2(.1);
plot(f, 'zebra')
axis equal off

%%
% Zebras, tigers, and bumblebees can be treated
% with adjustments of the colormap.
colormap([1 1 0; 0 0 0])
