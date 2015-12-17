%% Spike integral
% Nick Hale, October 2010

%%
% (Chebfun example quad/SpikeIntegral.m)
% [Tags: #quadrature, #spike]
%%
% We demonstrate the adaptive capabilities of Chebfun by integrating the
% "spike function"
f = @(x) sech(10*(x-0.2)).^2 + sech(100*(x-0.4)).^4 + ...
         sech(1000*(x-0.6)).^6 + sech(1000*(x-0.8)).^8;
%%
% (which appears as F21F in [1]) over $[0, 1]$.

%%
% The Chebfun representation is a very high degree polynomial,
% but this causes no difficulty.
ff = chebfun(f,[0 1])
LW = 'linewidth';
plot(ff,'b',LW,1.6), grid on
title('Spike function','FontSize',14)

%%
% Here is a confirmation that even the narrowest spike is well resolved:
semilogy(ff,'b','interval',[.795,.805],LW,1.6), grid on
title('Zoom, on semilogy axes','FontSize',14)

%%
% Now we compute the integral.  In order to estimate the time for this
% computation, we create the chebfun again without plotting it.
tic
ff = chebfun(f,[0 1]);
sum(ff)
toc

%%
% Now the degree of that polynomial was forced to be extraordinarily
% high in order to resolve the narrowest spike.  A much more
% compressed representation of $f$ can be attained by constructing the chebfun
% piecewise, using "splitting on".  As of December 2015,
% if this is done with default parameters, Chebfun fails to detect
% the narrowest spike:
ff = chebfun(f,[0 1],'splitting','on')
plot(ff,'b',LW,1.6), grid on
title('Unresolved spike function with splitting on','FontSize',14)

%%
% We can fix the problem by forcing Chebfun to sample at
% more points.  Note that the total number of parameters
% is 25 times less than with the global representation.
ff = chebfun(f,[0 1],'splitting','on','minSamples',100)
plot(ff,'b',LW,1.6), grid on
title('Resolved spike function with splitting on','FontSize',14)

%%
% If speed is all you care about, though, nothing has been gained over
% the first, global approach.  We compute the chebfun again
% and see that the integral is the same to full precision but
% the timing is worse:
tic
ff = chebfun(f,[0 1],'splitting','on','minSamples',100);
sum(ff)
toc

%% References
%
% 1. D. K. Kahaner, "Comparison of numerical quadrature formulas", in
%    J. R. Rice, ed., _Mathematical Software_, Academic Press, 1971, 229-259.
