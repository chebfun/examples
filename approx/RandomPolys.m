%% Roots of random polynomials
% Nick Trefethen, June 2014

%%
% (Chebfun example approx/RandomPolys.m)
% [Tags: #Das, #LEG2CHEB, #CHEB2LEG, #ROOTS]

%%
% Recently I heard a talk by Igor Pritsker of Oklahoma State
% University at which he discussed
% a theorem of Das in 1971 about the roots of random
% real polynomials [1].  This can be very nicely illustrated in Chebfun.

%%
% Das's result asserts that for a random polynomial with
% real coefficients, the fraction of roots that are real will be about
% $1/\sqrt 2 \approx 0.57735$.
% Specifically, this fraction appears (with probability $1$)
% asymptotically in the limit $n\to\infty$, where $n$ is the
% degree, assuming the random polynomials are defined as sums of Legendre
% polynomials with independent coefficients from the standard
% normal distribution.  The Legendre polynomials should be normalized
% by having 2-norm equal to $1$, not by taking the value $1$ at $x=1$.

%%
% Here for example is a random polynomial of degree 30:
rng('default');
n = 30;
cleg = randn(n+1,1);                           % Legendre coeffs
ccheb = leg2cheb(cleg,'normalization');        % Chebyshev coeffs
p = chebfun(ccheb,'coeffs');
LW = 'linewidth'; FS = 'fontsize'; MS = 'markersize';
plot(p,LW,1.6), axis([-1.1 1.1 -n n]), grid on
r = roots(p,'all');
rr = r(imag(r)==0);
hold on, plot(rr,p(rr),'.r',MS,16), hold off
ratio = length(rr)/n;
title(['fraction of real roots = ' num2str(ratio)],FS,12)

%%
% You might guess by looking at this picture that this was a
% polynomial of odd degree, but of course it is of even degree.
% So there must be a root outside $[-1,1]$, and we can confirm this
% by comparing
length(rr)
%%
% and 
length(roots(p))
%%
% Here is the rogue root:
rr(abs(rr)>1)

%%
% Now let's consider ten random polynomials of degree 500 and print
% the fraction of real roots for each:
n = 500;
data = [];
for k = 1:10
  cleg = randn(n+1,1);                         % Legendre coeffs
  ccheb = leg2cheb(cleg,'normalization');      % Chebyshev coeffs
  p = chebfun(ccheb,'coeffs');
  r = roots(p,'all');
  rr = r(imag(r)==0);
  ratio = length(rr)/n;
  data = [data ratio];
  disp(['fraction of real roots = ' num2str(ratio)])
end

%%
% The mean for the whole experiment is pretty close to $0.577$,
mean(data)

%% References
%
% 1. M. Das, Real zeros of a random sum of orthogonal polynomials,
%     _Proceedings of the American Mathematical Society_, 27
%     (1971), 147-153.
