%% Zeros of zeta(s) by analytic continuation
% Nick Trefethen and Mohsin Javed, July 2015

%%
% (Chebfun example complex/ZetaZeros.m)
% [Tags: #analytic continuation, #zeta]

%%
% The celebrated Riemann Hypothesis asserts that all the zeros of
% the zeta function (apart from those on the negative real axis) lie
% on the critical line $\hbox{Re} s = 1/2$ in the
% complex $s$-plane.  Computations of $\zeta(s)$
% and its zeros is a highly advanced subject and this Example certainly
% will not contribute to it!  However, we can show how easily certain
% kinds of analytic contination can be carried out in Chebfun.

%%
% A formula for $\zeta(s)$ that converges for $\hbox{Re}(s) >1$ is
% $$ \zeta(s) = \sum_{k=1}^\infty k^{-s}. $$
% For $\hbox{Re}(s) =4$, we can get approximately 16-digit
% precision with
% $$ \zeta(s) \approx \sum_{k=1}^{100000} k^{-s}. $$
% So here's our crude zeta function (note the summation in reverse
% order to minimize accumulation of rounding errors):
tic
zeta = @(s) sum((1e5:-1:1).^(-s));

%%
% For example, here are |zeta(4)| and the corresponding exact result:
zeta(4)
exact = pi^4/90

%%
% Let's work with a parameter $t \in [5, 50]$, and define
% $s = 4 + it$, so that 
% $s$ ranges over the complex interval $[4+5i, 4+50i]$.
s = chebfun(@(t) 4+1i*t,[5 50])

%%
% We can construct a chebfun (a complex function of the real parameter
% $t$) corresponding to the zeta function:
f = chebfun(@(t) zeta(s(t)),[5 50],'vectorize')

%%
% Here is the Chebfun ellipse of $f$ (see Chapter 8
% of [Trefethen 2013]) together with the numerically roots of $f$ in
% in the ellipse.  A black X is also marked to show the pole of
% the zeta function.
chebellipseplot(f), xlim([0 60]), axis equal, grid on
zeros_t = roots(f,'complex'); MS = 'markersize';
hold on, plot(zeros_t,'.r',MS,15)
plot(1i,3,'xk','markersize',8), hold off

%%
% Transplanted back to the $s$ variable, we see that the computed
% roots match the corresponding exact ones the 8 digits after the decimal point:
zeros_s = s(zeros_t);
zeros_exact = 0.5 + 1i*[14.1347251417 21.0220396388 25.0108575801 ...
30.4248761259 32.9350615877 37.5861781588 40.9187190121]';
ss = '%13.10f + %13.10fi   %13.10f + %13.10fi\n';
disp('            Chebfun                          Exact')
fprintf(ss,[real(zeros_s) imag(zeros_s) ...
            real(zeros_exact) imag(zeros_exact)].')

%%
% Along the critical line, here are the real and imaginary parts, with
% black dots showing the computed zeros.
t = chebfun('3.5i+t',[5 50]);
ft = f(t);
plot([imag(ft) real(ft)])
title('Real and imaginary parts of zeta(s) along critical line')
hold on, plot(real(zeros_t),imag(zeros_t-3.5i),'.k',MS,15)
grid on, hold off

%%
% Time taken by this example:
toc

%%
% Reference:
%
% 1. L. N. Trefethen, _Approximation Theory and Approximation
% Practice_, SIAM, 2013.
