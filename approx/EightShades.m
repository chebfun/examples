%% Eight Shades of Rational Approximation
% Mohsin Javed and Nick Trefethen, January 2016

%%
% (Chebfun example approx/EightShades.m)
% [Tags: #rational]

%% 1.  Introduction
% This example is a work in progress, not yet complete.

%%
% Our aim is to give a broad
% view of some practical methods of approximation of
% functions on an interval and of Chebfun's capabilities in
% these areas, not all of which are developed yet.
% As time goes by, this Example will probably evolve with
% new capabilities being described and new references given.

%%
% In a word, the "eight shades" come about as follows.
% We will discuss four types of approximation for nonperiodic
% functions (Chebyshev), and then their analogues for periodic
% functions (trigonometric).  Actually, each "four" is really
% a "four-and-a-half", since there is a least-squares option
% that inhabits the spectrum between the two extremes of
% interpolation (minimal number of data points) and projection
% (infinitely many data points).

%%
% And if you like it's not just eight or ten shades but
% sixteen or twenty! -- because for clarity, we begin by
% describing the better-known and better-developed polynomial
% special cases, i.e., type $(m,n)$ rational approxinations
% with $n=0$.

%% 2. Polynomial approximation
% If $f$ is a continuous function on $[-1,1]$, four
% interesting methods of approximation of $f$ by
% a degree $m$ polynomial are as follows.  These approximants
% can all be computed in Chebfun, and the mathematics
% is presented in _Approximation Theory and Approximation Practice_
% (ATAP).

%%
% _P1. Chebyshev interpolation_  (|chebfun| with |m+1| specified, _ATAP_ chap 4)
%
% _P2. Chebyshev projection_ (|chebfun| with |'trunc'| option, _ATAP_ chap 4)
% 
% _P3. Minimax approximation_ (|remez|, _ATAP_ chap 10)
%
% _P4. CF approximation_ (|cf|, _ATAP_ chap 20)

%%
% For example, here are degree 8 approximations of these kinds
% to $f(x) = \exp(-60(x-0.1)^2)$.
f = chebfun(@(x) exp(-50*(x-.1).^2),'trig'); m = 8;
FS = 'FontSize';
p1 = chebfun(f,m+1); subplot(2,2,1), yl = [-.5 1.2];
plot(f,'k',p1,'r'), ylim(yl), text(-.93,.9,'interpolation',FS,10)
p2 = chebfun(f,'trunc',m+1); subplot(2,2,2)
plot(f,'k',p2,'r'), ylim(yl), text(-.93,.9,'projection',FS,10)
p3 = remez(f,m); subplot(2,2,3)
plot(f,'k',p3,'r'), ylim(yl), text(-.93,.9,'minimax',FS,10)
p4 = cf(f,m); subplot(2,2,4)
plot(f,'k',p4,'r'), ylim(yl), text(-.93,.9,'CF',FS,10)

%%
% These curves show some properties that are typical of
% such approximations.  One is that the
% differences between them are not very great.
% Another (for smooth functions $f$, at least) is that 
% the minimax and CF approximations, though mathematically distinct,
% are for practical purposes indisinguishable.  We can quantify
% this effect for the present example by measuring the maximal
% difference between the two:
CFerror = norm(p3-p4,inf)

%%
% Methods P1 and P2 represent two ends of a spectrum.  In between,
% there is a method we could label P1.5:

%%
% _P1.5. Chebyshev least-squares_  (|ratinterp|, _ATAP_ chap 26)

%%
% The idea here is to determine a polynomial $p$ of the specified
% degree $m$ that is the least-squares approximation to $f$
% on the $K$-point Chebyshev grid, where $K$ satisfies
% $m+1 \le K < \infty$.  For $K=m+1$, this is the same as
% Chebyshev interpolation, and in the limit $K \to \infty$
% the discrete least-squares problem becoming a continuous
% least-squares problem with the Chebyshev weight, i.e., Chebyshev
% projection.  Chebfun has no special code for computing
% Chebyshev least-squares approximation polynomials;
% the code |ratinterp| and the _ATAP_ chapter cited above
% both apply more generally to the rational case.

%% 3. Trigonometric polynomial approximation
% The four methods of polynomial approximation have trigonometric
% analogues for periodic functions.  Our favorite starting reference on
% this material is [2].  At present, Chebfun has
% a |trigremez| command for trigonometric minimax approximation, but
% not yet a |trigcf| command.  A |triginterp| command for the
% least-squares case is under development but not yet in the
% development or master branches of Chebfun.

%%
% _TP1. Trigonometric interpolation_ (|chebfun| with |'trig'| specified, [2])
%
% _TP2. Trigonometric projection_ (|chebfun| with |'trunc'| and |'trig'|, [2])
% 
% _TP3. Minimax trigonometric approximation_ (|trigremez|)
%
% _TP4. Fourier-CF approximation_ (|trigcf|, not yet available)

%%
% Again, TP1 and TP2 represent two ends of a spectrum: 

%%
% _TP1.5. Trigonometric least-squares_ (|triginterp|, in a branch)

clf
t1 = chebfun(f,m+1,'trig'); subplot(2,2,1)
plot(f,'k',t1,'b'), ylim(yl), text(-.93,.9,'interpolation',FS,10)
t2 = chebfun(f,'trunc',m+1,'trig'); subplot(2,2,2)
plot(f,'k',t2,'b'), ylim(yl), text(-.93,.9,'projection',FS,10)
t3 = trigremez(f,m/2); subplot(2,2,3)
plot(f,'k',t3,'b'), ylim(yl), text(-.93,.9,'minimax',FS,10)
subplot(2,2,4)
text(-.93,.9,'CF',FS,10)
text(-.5,.2,'(not yet available)',FS,10), axis([-1 1 yl])
set(gca,'xtick',[],'ytick',[])

%% 4. Rational approximation

%%
% Discussion to be added here.

%%
% _R1. Rational interpolation_ (|ratinterp|, chap 27 of _ATAP_)
%
% _R2. Chebyshev-Pade approximation_ (|chebpade|)
% 
% _R3. Minimax rational approximation_ (|remez|, chap 24 of _ATAP_)
%
% _R4. CF rational approximation_  (|cf|, chap 20 of _ATAP_)

%%
% and

%%
% _R1.5. Rational least-squares_ (|ratinterp|, chap 27 of _ATAP_)

clf
m = 3; n = 3;
[p,q] = ratinterp(f,m,n); r1 = p./q; subplot(2,2,1), yl = [-.5 1.2];
plot(f,'k',r1,'r'), ylim(yl), text(-.93,.9,'interpolation',FS,10)
[p,q] = chebpade(f,m,n); r2 = p./q; subplot(2,2,2)
plot(f,'k',r2,'r'), ylim(yl), text(-.93,.9,'projection',FS,10)
[p,q] = remez(f,m,n); r3 = p./q; subplot(2,2,3)
plot(f,'k',r3,'r'), ylim(yl), text(-.93,.9,'minimax',FS,10)
[p,q] = cf(f,m,n); r4 = p./q; subplot(2,2,4)
plot(f,'k',r4,'r'), ylim(yl), text(-.93,.9,'CF',FS,10)

%% 5. Trigonometric rational approximation

%%
% Discussion to be added here.

%%
% _TR1. Trigonometric rational interpolation_ (|triginterp|, in a branch)
%
% _TR2. Fourier-Pade approximation_ (|trigpade|, in a branch)
% 
% _TR3. Minimax trigonometric rational approx_ (|trigremez|, not yet available)
%
% _TR4. Fourier-CF rational approximation_  (|trigcf|, not yet available)

%%
% and

%%
% _TR1.5. Trigonometric rational least-squares_  (|triginterp|, in a branch)

clf
subplot(2,2,1)
text(-.93,.9,'interpolation',FS,10)
text(-.5,.2,'(not yet available)',FS,10), axis([-1 1 yl])
set(gca,'xtick',[],'ytick',[])
subplot(2,2,2)
text(-.93,.9,'projection',FS,10)
text(-.5,.2,'(not yet available)',FS,10), axis([-1 1 yl])
set(gca,'xtick',[],'ytick',[])
subplot(2,2,3)
text(-.93,.9,'minimax',FS,10)
text(-.5,.2,'(not yet available)',FS,10), axis([-1 1 yl])
set(gca,'xtick',[],'ytick',[])
subplot(2,2,4)
text(-.93,.9,'CF',FS,10)
text(-.5,.2,'(not yet available)',FS,10), axis([-1 1 yl])
set(gca,'xtick',[],'ytick',[])

%% 6. References
%
% [1] L. N. Trefethen, _Approximation Theory and
% Approximation Practice_, SIAM, 2013.
%
% [2] G. B. Wright, M. Javed, H. Montanelli and L. N. Trefethen,
% Extension of Chebfun to periodic functions, _SIAM J. Sci. Comp._,
% 2016.
