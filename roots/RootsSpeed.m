%% Speed and Accuracy of Chebfun `roots`
% Jared Aurentz and Nick Trefethen, 21st May 2014

%%
% (Chebfun example roots/RootsSpeed.m)
% [Tags: #colleaguematrix, #ROOTS]

%% 1. What Chebfun does
% Let us do a test of the speed and accuracy of the
% Chebfun algorithm for computing roots.  We'll pick
% a function whose roots we know:
f = chebfun(@(x) exp(x).*sin(1000*pi*x));

%%
% Here is the length of the function, approximately $1000\pi$.
n = length(f)

%%
% This function has 2001 roots linearly spaced from
% $-1$ to $1$:
exact = linspace(-1,1,2001)';

%%
% Let's compute them with Chebfun.  We'll do it
% twice just to make sure the timing estimate is
% realistic.
r = roots(f);
tic, r = roots(f); 

%%
% Here is the time elapsed:
toc

%%
% Here is the maximum error:
norm(r-exact,inf)

%%
% This number is twice machine epsilon, which looks
% very good.  (However, you couldn't ask for a better
% conditioned problem than this one, since the derivative
% of $f$ at each root is large.)

%%
% For fun we plot $f$ and its roots over a short interval:
d = [-0.0105,0.0105]; 
plot(f,'interval',d,'linewidth',1.6)
axis([d -1.5 1.5]), grid on
hold on, plot(r,f(r),'.r','markersize',18), hold off


%% 2. What Chebfun might do
% "Classically", using the MATLAB `roots` command,
% it takes $O(n^3)$ operations to compute the
% eigenvalues of a companion matrix, which is the
% method that MATLAB has used since the 1970s for
% finding roots of a polynomial in the monomial
% basis.  We can illustrate that this computation is
% slow, if not really that the complexity is cubic,
% by the following experiment.
for ntest = [250 500 1000]
    c = randn(ntest,1);
    tic, roots(c); toc
end

%%
% It is clear from this experiment that calling
% Matlab `roots` for a polynomial of degree as large
% as our chebfun `f` would be very slow.

%%
% However, $O(n^2)$ algorithms for this
% problem have been available for quite a while,
% though they are not built into MATLAB.
% A key person in this effort over the years has
% been Dario Bini of the
% University of Pisa.  See for example
% [Bini et al. 2010].  This group also offers Fortran
% software.  Another notable contribution is
% [Chandrasekaran et al. 2007].

%%
% The first author of this example, in collaboration
% with David Watkins and others, has been developing alternative
% $O(n^2)$ algorithms for the companion matrix eigenvalue
% problem [Aurentz et al. 2014].

%%
% All this is for the companion matrix eigenvalue problem,
% which corresponds to polynomial rootfinding in the monomial
% basis, a problem that makes sense when your roots are on or
% near the unit circle.  What about the Chebfun context of
% roots on or near $[-1,1]$?  Here the analogous matrix structure
% is a so-called _colleague matrix_, dating to Specht and Good
% around 1960; see Chapter 18 of [Trefethen 2013].  What
% can be done in this case?

%%
% First of all we note what Chebfun currently does: following
% an idea of Boyd [Boyd 2002], it subdivides
% the interval when necessary.  This is how a superficially
% $O(n^3)$ algorithm is brought down to $O(n^2)$, enabling the
% good performance seen above.  Note that intervals are
% different from circles: if you split an interval in half, you
% get two intervals, but if you split a circle in half, you don't
% get two circles.  Therefore this recursion trick is not
% available in the monomial case.

%%
% But what about an $O(n^2)$ _linear algebra_ solution to the
% problem, rather than relying on splitting of intervals?
% Here too there has been progress [Eidelman et al. 2008], and
% we hope to have work of our own to report before long.
% An interesting project for the future will be to see whether
% Chebfun rootfinding can be improved by the use of $O(n^2)$ linear algebra
% algorithms rather than interval subdivision,
% while holding to Chebfun's principle of doing everything
% in MATLAB without relying on Mex files to link to Fortran or C.

%% References
%
% [Aurentz et al. 2014]
% J. Aurentz, T. Mach, R. Vandebril and D. S. Watkins,
% to appear.
%
% [Bini et al. 2010]
% D. A. Bini, P. Boito, Y. Eidelman, L. Gemignani, and
% I. Gohberg, A fast implicit QR eigenvalue algorithm
% for companion matrices, _Linear Algebra and its
% Applications_, 432 (2010), 2006-2031.
%
% [Boyd 2002] J. P. Boyd, Computing zeros on a real interval
% through Chebyshev expansion and polynomial rootfinding,
% _SIAM Journal on Numerical Analysis_, 40 (2002), 1666-1682.
%
% [Chandrasekaran et al. 2007]
% S. Chandrasekaran, M. Gu, J. Xia and J. Zhu,
% A fast QR algorithm for companion matrices,
% _Operator Theory: Advances and Applications_, 179
% (2007), 111-143.
%
% [Eidelman et al. 2008] Y. Eidelman, L. Gemignani, and
% I. Gohberg, Efficient eigenvalue computation for
% quasiseparable Hermitian matrices under low
% rank perturbations, _Numerical Algorithms_,
% 47 (2008), 253-273.
%
% [Trefethen 2013]
% L. N. Trefethen, _Approximation Theory and Approximation
% Practice_, SIAM, 2013.

