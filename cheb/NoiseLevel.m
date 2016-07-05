%% What's the noise level of a chebfun?
% Nick Trefethen, July 2016

%%
% (Chebfun example cheb/NoiseLevel.m)

%% 1. The problem of estimating noise level
% A familar image for any Chebfun user is the ``noise plateau''
% that results from rounding errors in a Chebyshev series.
% To see this, it's enough to construct a chebfun with 
% the `doublelength` flag, like this:
f = chebfun(@(x) cos(100*exp(x)),'doublelength');
plotcoeffs(f,'.'), xlim([0 450])

%%
% The question considered in this memo is, given an
% already constructed chebun, can we estimate its
% noise level?  For example suppose we make a
% chebfun from coefficients 0 through 280 of the
% series just plotted:
c = chebcoeffs(f);
f2 = chebfun(c(1:281),'coeffs');
plotcoeffs(f2,'.'), xlim([0 450])

%% 
% You and I know that the last 60 or so dots are just noise.
% However, Chebfun's `simplify` command at present (July 5, 2016)
% fails to detect this properly:
plotcoeffs(simplify(f2),'.'), xlim([0 450])

%%
% In Section 2 I'll describe the method currently used by `simplify`,
% and in Section 3 I'll propose an alternative method.  Section 4
% will put all the in the wider context of Chebfun computations.

%% 2. Estimating noise via FFT and IFFT
% OK, suppose we're given a Chebyshev series
% like that of f2 above.  Here is
% one way to estimate noise: pad it with zeros, then take a Chebyshev
% transform followed by an inverse Chebyshev transform.  Rounding errors
% in the transforms will turn your zeros into something nonzero,
% and that's your estimated noise level.

%%
% Let's try it.  We pad with 150 zeros, transform, and transform back.
c2 = chebcoeffs(f2); c2 = [c2; zeros(150,1)];
v2 = chebtech2.coeffs2vals(c2); c3 = chebtech2.vals2coeffs(v2);
MS = 'markersize';
semilogy(0:length(c3)-1,abs(c3),'or',MS,3), hold on
plotcoeffs(f2,'.'), xlim([0 450]), hold off

%%
% This is not so good!  The red dots lie about two orders
% of magnitude too low.  This is why `simplify` has gotten
% the wrong answer.  It has estimated the noise
% level in this fashion to be around $10^{-17}$, relative to
% which the coefficients of size $10^{-15}$ appear to be
% genuine.  That's why it hasn't discarded those coefficients.

%% 3. Estimating noise with x-perturbations
% What went wrong with the idea of transforming forward and then
% back again?  I think the trouble is that this way of thinking
% takes into account ``vertical rounding errors but not horizontal
% ones''.  By horizontal rounding errors, I mean those associated
% with errors in the argument $x$ to a function like $f(x)$.  
% Because of the $\cos(100x)$ factor, our function has derivatives
% of size $O(100)$, 

%% 4. The wider Chebfun context
% During 2013-2015, each chebfun
% contained an estimate of this noise level as a guide
% to further computations.  In the end we abandoned this
% idea because it lacked a precise definition and too
% often led to non-optimal decisions.  Instead, in 2015
% we introduced series chopping by the method encoded
% in `standardChop.m`, described in [1], which 

%% 5. References
%
% 1. J. L. Aurentz and L. N. Trefethen, Chopping a Chebyshev
% series, to appear in _ACM Transactions on Mathematical Software_.
