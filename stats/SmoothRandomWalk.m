%% Smooth random walk
% Nick Trefethen, February 2017

%%
% (Chebfun example/stats/SmoothRandomWalk.m )
% [Tags: #randnfun]

%%  
% By integrating coin flips in one or more dimensions,
% we get a random walk, which becomes Brownian motion
% in the limit of infinitely many infinitely small steps.
% Chebfun's |randnfun| command enables us to explore
% a smooth continuous analogue of this process.

%%
% Let's work in 2D, using a complex variable for convenience.
% Here we plot the indefinite integral of a complex random function scaled
% by $(dx)^{-1/2}$.  Red dots mark the initial and end points.
LW = 'linewidth'; FS = 'fontsize'; MS = 'markersize';
dx = 0.1;
rng(3)
f = (randnfun(dx) + 1i*randnfun(dx))/sqrt(dx);
g = cumsum(f);
plot(g,'k',LW,1), grid on, hold on
plot(g([-1 1]),'.r',MS,8), hold off
axis(1.5*[-1 1 -1 1]), axis square
title(['dx = ' num2str(dx)],FS,14)

%%
% We divide the characteristic length defining 
% the random function by 4 three times.
% The limit of Brownian motion is being approached,
% as no doubt could be proved.

for k = 1:3
  dx = dx/4;
  f = (randnfun(dx) + 1i*randnfun(dx))/sqrt(dx);
  g = cumsum(f);
  plot(g,'k',LW,1), grid on, hold on
  plot(g([-1 1]),'.r',MS,8), hold off
  axis(1.5*[-1 1 -1 1]), axis square
  title(['dx = ' num2str(dx)],FS,14), snapnow
end

%%
% Here is a zoom of the final image:
axis([-1.4 .2 -.8 .8]), axis square, axis off, title(' ')
