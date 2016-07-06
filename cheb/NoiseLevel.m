%% What's the noise level of a chebfun?
% Nick Trefethen, 6 July 2016

%%
% (Chebfun example cheb/NoiseLevel.m)

%% 1. The problem of estimating noise level
% A familar image for any Chebfun user is the ``noise plateau''
% that results from rounding errors in a Chebyshev series.
% To see this, it's enough to construct a chebfun with 
% the `doublelength` flag, like this:
%%
% <latex> \vskip -2em </latex>
f = chebfun(@(x) cos(100*exp(x)),'doublelength');
MS = 'markersize';
plotcoeffs(f,'.b',MS,5), axis([0 450 1e-18 10])
set(gca,'ytick',10.^(-15:5:0))
%%
% <latex> \vskip .1em </latex>

%%
% The question to be considered here is, given an
% already constructed chebfun, can we estimate its
% noise level?  For example, suppose we make a
% chebfun from coefficients 0 through 280 of the
% series just plotted:
%%
% <latex> \vskip -2em </latex>
c = chebcoeffs(f);
f2 = chebfun(c(1:281),'coeffs');
plotcoeffs(f2,'.b',MS,5), axis([0 450 1e-18 10])
set(gca,'ytick',10.^(-15:5:0))
%%
% <latex> \vskip .1em </latex>

%% 
% You and I know that the last 60 or so dots are just noise.
% However, Chebfun's `simplify` command at present (July 6, 2016)
% fails to detect this properly.  We call simplify, yet
% nothing changes:
%%
% <latex> \vskip -2em </latex>
plotcoeffs(simplify(f2),'.b',MS,5), axis([0 450 1e-18 10])
set(gca,'ytick',10.^(-15:5:0))
%%
% <latex> \vskip .1em </latex>

%%
% Why hasn't Chebfun gotten this right?

%%
% In Section 2 I'll describe the method currently used by `simplify`,
% and then in Section 3 I'll propose an alternative.  Section 4
% gives more examples, and Section 5 puts
% the discussion in the wider context of Chebfun computations.

%% 2. Estimating noise via FFT and IFFT
% OK, suppose we're given a Chebyshev series
% like that of `f2` above.  Here is the method used by Chebfun `simplify`
% to estimate the noise level:
% pad the series with zeros, then take a Chebyshev
% transform followed by an inverse Chebyshev transform (essentially,
% FFT and inverse FFT).  Rounding errors
% in the transforms will turn the zeros into nonzeros,
% and that's the estimated noise level.

%%
% Let's try it.  Here is a short code called `extend1` to pad
% a series with a bunch of zeros, transform, and transform
% back.  
%%
% <latex> \vskip -2em </latex>
type extend1

%%
% Here we apply `extend1` to the data above.
%%
% <latex> \vskip -2em </latex>
c2 = chebcoeffs(f2);
c3 = extend1(c2);
semilogy(0:length(c3)-1,abs(c3),'or',MS,1), hold on
plotcoeffs(f2,'.b',MS,5), axis([0 450 1e-18 10]), hold off
set(gca,'ytick',10.^(-15:5:0))
%%
% <latex> \vskip .1em </latex>

%%
% This is not good!  The red dots lie two orders
% of magnitude too low.  This is why `simplify` has gotten
% the wrong answer.  It has estimated the noise
% level for this function to be around $10^{-17}$, and based on that estimate,
% it has judged that the coefficients of size $10^{-15}$ are
% genuine.  That's why it hasn't discarded those coefficients.
% (The ``it'' in this discussion is actually the code `standardChop`.)

%% 3. Including x-perturbations
% What went wrong with the idea of transforming forward and then
% back again?  I think the trouble is that this way of thinking
% takes into consideration ``vertical'' rounding errors but not ``horizontal''
% ones.  By horizontal rounding errors, I mean those associated
% with perturbations on the order of machine epsilon
% in the argument $x$.  
% For our particular function $f(x) = \cos(100e^x)$, the derivatives
% will be of size $O(100)$, so these errors will be not
% of order machine precision but $O(100)$ times as large.

%%
% You might think that tracking this would be a complicated matter
% involving derivative estimates, but I think something simpler can
% do the trick: after transforming from coefficients to data values,
% perturb the data values at random ``to the left and right''.
% What I mean by this is, rather than take a data value straight, blend
% in a very small component of the value to its left or right, with
% the blending coefficients determined in the natural way based on
% the ratio of machine epsilon to the local grid spacing.
% Here is a code that does this, not very compact I fear, subject
% to improvement.
%%
% <latex> \vskip -2em </latex>
type extend2

%%
% If we apply `extend2` instead of `extend1` to the data above, the 
% result is much better.  In all the plots, red dots come from
% `extend1` and green dots from `extend2`.
%%
% <latex> \vskip -2em </latex>
c2 = chebcoeffs(f2);
c4 = extend2(c2);
semilogy(0:length(c4)-1,abs(c4),'og',MS,1), hold on
plotcoeffs(f2,'.b',MS,5), axis([0 450 1e-18 10]), hold off
set(gca,'ytick',10.^(-15:5:0))
%%
% <latex> \vskip .1em </latex>

%% 4. Some more examples
% Let's try some more functions.  In each case, the result from
% `extend1` is on the left with the red dots, and the result from
% `extend2` is on the right with the green dots.

%%
% <latex> \vskip -2em </latex>
f = chebfun(@(x) cos(1000*exp(x)),'doublelength');
c = chebcoeffs(f); c = c(1:round(.65*length(f)));
subplot(1,2,1), d = extend1(c);
semilogy(0:length(d)-1,abs(d),'or',MS,1), hold on
semilogy(abs(c),'.b',MS,5), hold off
subplot(1,2,2), d = extend2(c);
semilogy(0:length(d)-1,abs(d),'og',MS,1), hold on
semilogy(abs(c),'.b',MS,5), hold off
%%
% <latex> \vskip .1em </latex>

%%
% The next example shows an effect that's surprising at
% first -- why are there green dots to the left, where there
% are no red dots?  What's going on is that this is an
% even function, with all its odd-order coefficients zero and
% off the bottom of the scale;
% `extend1` maintains that symmetry whereas `extend2` does
% not.  It would be easy enough to make `extend2` preserve
% symmetry if that were considered desirable.

%%
% <latex> \vskip -2em </latex>
f = chebfun(@(x) 1./(1+25*x.^2),'doublelength');
c = chebcoeffs(f); c = c(1:round(.65*length(f)));
subplot(1,2,1), d = extend1(c);
semilogy(0:length(d)-1,abs(d),'or',MS,1), hold on
semilogy(abs(c),'.b',MS,5), hold off, ylim([1e-20 1])
subplot(1,2,2), d = extend2(c);
semilogy(0:length(d)-1,abs(d),'og',MS,1), hold on
semilogy(abs(c),'.b',MS,5), hold off, ylim([1e-20 1])
%%
% <latex> \vskip .1em </latex>

%%
% <latex> \vskip -2em </latex>
f = chebfun(@(x) 1./(1+25*x.^2));
c = chebcoeffs(f); c = c(1:round(.5*length(f)));
subplot(1,2,1), d = extend1(c);
semilogy(0:length(d)-1,abs(d),'or',MS,1), hold on
semilogy(abs(c),'.b',MS,5), hold off
subplot(1,2,2), d = extend2(c);
semilogy(0:length(d)-1,abs(d),'og',MS,1), hold on
semilogy(abs(c),'.b',MS,5), hold off
%%
% <latex> \vskip .1em </latex>

%%
% <latex> \vskip -2em </latex>
f = chebfun(@(x) 1./(1+2500*x.^2),'doublelength');
c = chebcoeffs(f); c = c(1:round(.65*length(f)));
subplot(1,2,1), d = extend1(c);
semilogy(0:length(d)-1,abs(d),'or',MS,1), hold on
semilogy(abs(c),'.b',MS,5), hold off, ylim([1e-22 1])
subplot(1,2,2), d = extend2(c);
semilogy(0:length(d)-1,abs(d),'og',MS,1), hold on
semilogy(abs(c),'.b',MS,5), hold off, ylim([1e-22 1])
%%
% <latex> \vskip .1em </latex>

%%
% <latex> \vskip -2em </latex>
f = chebfun(@(x) log(1.0001-x),'doublelength');
c = chebcoeffs(f); c = c(1:round(.65*length(f)));
subplot(1,2,1), d = extend1(c);
semilogy(0:length(d)-1,abs(d),'or',MS,1), hold on
semilogy(abs(c),'.b',MS,5), hold off, ylim([1e-20,1e2])
subplot(1,2,2), d = extend2(c);
semilogy(0:length(d)-1,abs(d),'og',MS,1), hold on
semilogy(abs(c),'.b',MS,5), hold off, ylim([1e-20,1e2])
%%
% <latex> \vskip .1em </latex>

%%
% <latex> \vskip -2em </latex>
f = chebfun(@(x) exp(-100./x.^2),'doublelength');
c = chebcoeffs(f); c = c(1:round(.65*length(f)));
subplot(1,2,1), d = extend1(c);
semilogy(0:length(d)-1,abs(d),'or',MS,1), hold on
semilogy(abs(c),'.b',MS,5), hold off, ylim([1e-65 1e-40])
subplot(1,2,2), d = extend2(c);
semilogy(0:length(d)-1,abs(d),'og',MS,1), hold on
semilogy(abs(c),'.b',MS,5), hold off, ylim([1e-65 1e-40])
%%
% <latex> \vskip .1em </latex>

%%
% <latex> \vskip -2em </latex>
f = chebfun(@(x) exp(-100./(x+1).^2),'doublelength');
c = chebcoeffs(f); c = c(1:round(.65*length(f)));
subplot(1,2,1), d = extend1(c);
semilogy(0:length(d)-1,abs(d),'or',MS,1), hold on
semilogy(abs(c),'.b',MS,5), hold off
subplot(1,2,2), d = extend2(c);
semilogy(0:length(d)-1,abs(d),'og',MS,1), hold on
semilogy(abs(c),'.b',MS,5), hold off
%%
% <latex> \vskip .1em </latex>

%% 5. The wider Chebfun context
% During 2013-2015, each chebfun
% contained an estimate of this noise level as a guide
% to further computations.  In the end we abandoned this
% idea because it lacked a precise definition and too
% often led to non-optimal decisions.  Instead, in 2015
% we introduced series chopping by the method encoded
% in `standardChop.m`, described in [1], which 
% doesn't use noise estimates anywhere _except_ in
% the command simplify.  So the whole subject of this
% memo lies not at the heart of Chebfun (namely construction)
% but only in the more specialized area of simplification.

%% 6. Does it work?
% Does the new method eliminate anomalies?  Well, let's see.
% At the beginning we constructed a function and simplified
% it with no effect.  Here's that image again:
%%
% <latex> \vskip -2em </latex>
f = chebfun(@(x) cos(100*exp(x)),'doublelength');
c = chebcoeffs(f);
f2 = chebfun(c(1:281),'coeffs');
f3 = simplify(f2);
clf, plotcoeffs(f3,'.b',MS,5), axis([0 450 1e-18 10])
set(gca,'ytick',10.^(-15:5:0))
%%
% <latex> \vskip .1em </latex>

%%
% Now we do it again, but this time with a call to `simplify2`,
% a modified version of `simplify` based on `extend2`.
%%
% <latex> \vskip -2em </latex>
f2 = chebfun(c(1:281),'coeffs');
f4 = simplify2(f2);
plotcoeffs(f4,'.b',MS,5), axis([0 450 1e-18 10])
set(gca,'ytick',10.^(-15:5:0))
%%
% <latex> \vskip .1em </latex>

%%
% Another place where simplify has caused trouble is in Chebfun3.
% For example, the lengths 4097 in this construction are suspicious:
%%
% <latex> \vskip -2em </latex>
f = chebfun3(@(x,y,z) cos(1000*pi*(x+y+z)))
%%
% <latex> \vskip .1em </latex>

%%
% Sure enough, `plotcoeffs` confirms that the series have
% not been chopped properly.
%%
% <latex> \vskip -2em </latex>
plotcoeffs(f)
%%
% <latex> \vskip .1em </latex>

%%
% I haven't yet tested whether `simplify2` will fix this problem too.


%% 7. Reference
%
% 1. J. L. Aurentz and L. N. Trefethen, Chopping a Chebyshev
% series, to appear in _ACM Transactions on Mathematical Software_.
