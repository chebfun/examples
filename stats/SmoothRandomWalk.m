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
lw = 2.5; fs = 32; ms = 40;
dx = 0.1;
rng(5)
f = randnfun(dx,'norm') + 1i*randnfun(dx,'norm');
g = cumsum(f);
plot(g,'k',LW,lw), grid on, hold on
plot(g([-1 1]),'.r',MS,ms), hold off
axis(1.9*[-1 1 -1 1]), axis square
title(['dx = ' num2str(dx)],FS,fs)
set(gca,'xtick',-2:2,'ytick',-2:2)

%%
% We divide the characteristic length defining 
% the random function by 4 three times.
% The limit of Brownian motion is being approached,
% as no doubt could be proved.

for k = 1:3
  dx = dx/4;
  f = randnfun(dx,'norm') + 1i*randnfun(dx,'norm');
  g = cumsum(f);
  plot(g,'k',LW,lw), grid on, hold on
  plot(g([-1 1]),'.r',MS,ms), hold off
  axis(1.9*[-1 1 -1 1]), axis square
  title(['dx = ' num2str(dx)],FS,fs)
  set(gca,'xtick',-2:2,'ytick',-2:2), snapnow
end

%%
% Here is a zoom of the final image:
axis([-1.6 0.6 -1.1 1.1]), axis square, axis off, title(' ')
snapnow
