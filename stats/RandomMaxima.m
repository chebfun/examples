%% Maxima of a random function
% Nick Trefethen, February 2017

%%
% (Chebfun example/stats/RandomMaxima.m )
% [Tags: #randnfun]

%%  
% How many maxima does a random function have?
% We can look at one answer to this question based on
% the notion of a random function implemented in the
% new Chebfun command |randnfun|.  Here for example are
% the local maxima of a random function on the interval $[0,20]$.
tic, rng(0), dx = 1;
f = randnfun(dx,[0,20]);
[val,pos] = max(f,'local');
MS = 'markersize'; LW = 'linewidth'; FS = 'fontsize';
plot(f,'k',LW,1.6), grid on, hold on
plot(pos,val,'.r',MS,16), hold off
title([int2str(length(val)) ' maxima'],FS,16)

%%
% The random function in question is obtained as the restriction of a
% sum of a finite Fourier series with independent normally distributed
% random coefficients; the maximum wave number is about $2\pi /dx$, i.e.,
% $2\pi$ in this case with $dx = 1$.
% Here is an analogous picture on an interval of length 40.
f = randnfun(dx,[0,40]);
[val,pos] = max(f,'local');
plot(f,'k',LW,1.6), grid on, hold on
plot(pos,val,'.r',MS,16), hold off
title([int2str(length(val)) ' maxima'],FS,16)

%%
% Let's explore the number of maxima as a function of the length
% of the interval:
Lvec = 2.^(0:10);
nmaxvec = [];
for L = Lvec
  f = randnfun(dx,[0 L]);
  nmaxvec = [nmaxvec length(max(f,'local'))];
end 
loglog(Lvec,Lvec,'-r'), hold on
loglog(Lvec,nmaxvec,'.',MS,20), grid on, hold off
axis([.8 1300 .8 1300]); 
xlabel('length of interval',FS,16)
ylabel('no. of maxima',FS,16)

%%
% It would seem that the expected number of maxima is asymptotic to $L$.
% Very likely there is a literature on this question.

%%
% Note that the first two dots seem higher than the curve.
% This is because when Chebfun finds "local extrema", it includes
% extrema at the endpoints even though these (with probability 1) will
% not be points of zero derivative.  We could improve the experiment
% to correct for this, though this is not too important since the
% effect diminishes proportionally when $L$ is large.  

%%
% Execution time for this example:
format short, time_in_seconds = toc
