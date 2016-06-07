%% Low-rank compression of square and round pegs
% Nick Trefethen, June 2016

%%
% (Chebfun example approx2/Pegs.m)

%%
% A recent example called "Low-rank approximation and
% alignment with axes" explored the phenomenon that
% the alignment of a function $f(x,y)$ with the $x$ and $y$
% axes may greatly affect the rank of its representation
% in Chebfun2.  Here we look at three more examples from
% `cheb.gallery2`.

%%
% Here's the "tilted peg" example.
[f,fa] = cheb.gallery2('tiltedpeg'); fa

%%
% A plot shows this this is a somewhat smoothed version of
% the characteristic function of a tilted square:
levels = .1:.2:.9
contourf(f,levels), axis equal
set(gca,'xtick',-1:1,'ytick',-1:1)
FS = 'fontsize';
text(-.9,.8,['rank ' int2str(rank(f))],FS,18)

%%
% If the peg is aligned with the axes it is a separable
% function, hence of rank 1:
[f,fa] = cheb.gallery2('squarepeg'); fa
contourf(f,levels), axis equal
set(gca,'xtick',-1:1,'ytick',-1:1)
text(-.9,.8,['rank ' int2str(rank(f))],FS,18)

%%
% A round peg has an in-between rank:
[f,fa] = cheb.gallery2('roundpeg'); fa
contourf(f,levels), axis equal
set(gca,'xtick',-1:1,'ytick',-1:1)
text(-.9,.8,['rank ' int2str(rank(f))],FS,18)

%%
% There are analogous examples in the Diskfun gallery
% collection, to be released before long.  There, the
% round peg has rank 1 whereas the square peg does not.
% Tilting does not matter.
