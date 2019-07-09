%% 2D zero set example of Dmitry Belyaev
% Nick Trefethen, July 2019

%%
% (Chebfun example approx2/Belyaev.m)

%%
% Dmitry Belyaev at Oxford is an expert on zero sets of
% functions composed form random plane waves and related
% problems.  Here is an example he has looked at:
tic
LW = 'linewidth'; XT = 'xtick'; YT = 'ytick'; FS = 'fontsize';
rng(1); a = randn(1,4) + 1i*rand(1,4);
cheb.xy
wave = @(k) real(a(1)*exp(i*pi*(k*x-y)) + a(2)*exp(i*pi*(k*x+y)) ...
               + a(3)*exp(i*pi*(k*y-x)) + a(4)*exp(i*pi*(k*y+x)));
r = roots(wave(8));
plot(r, LW, 2)
axis([-1 1 -1 1]), axis square, set(gca,XT,[],YT,[])
title(['number of components: ' int2str(size(r,2))],FS,10)

%%
% The Chebfun2 |roots| command has picked out the distinct
% components of the zero set in the unit square.  These computations
% are delicate, however, and Chebfun does not always get this number
% right (nor does it always get the curves accurately).
% Here we set $k=16$:
r = roots(wave(16));
plot(r, LW, 1.2)
axis([-1 1 -1 1]), axis square, set(gca,XT,[],YT,[])
title(['number of components: ' int2str(size(r,2))],FS,10)

%%
% And here we set $k=32$:
r = roots(wave(32));
plot(r, LW, .7)
axis([-1 1 -1 1]), axis square, set(gca,XT,[],YT,[])
title(['number of components: ' int2str(size(r,2))],FS,10)

%%
% Total time for this example:
toc
