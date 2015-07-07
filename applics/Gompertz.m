%% Exponential, logistic, and Gompertz growth
% Toby Driscoll, June 15, 2015

%%
% (Chebfun example applics/Gompertz.m)
% [Tags: #ODE, #IVP]

%%
% If the per-capita growth rate of a population is held constant, exponential
% growth of the population results. 
exponentialPCRate = @(P) 0.5;
N = chebop( @(t,P) diff(P) - P.*exponentialPCRate(P), [0 25] );
N.lbc = 0.2;
exponentialResult = N\0;
plot( exponentialResult )
xlabel('t'), ylabel('P(t)'), title('Exponential growth')

%% 
% The result is unbounded growth, which is not biologically realistic. A
% common way to remedy this defect is the *logistic model*. In this model,
% the per capita growth rate decreases linearly to zero as the population
% $P$ approaches a fixed value, known as the *carrying capacity*. For
% example, a carrying capacity of $P=6$ is imposed through
logisticPCRate = @(P) 0.5*(6-P)/5.8;

%%
% Here is the resulting growth.
N = chebop( @(t,P) diff(P) - P.*logisticPCRate(P), [0 25] );
N.lbc = 0.2;
logisticResult = N\0;
plot( [exponentialResult, logisticResult] )
ylim([0 8])
xlabel('t'), ylabel('P(t)'), legend('exp','logistic')

%%
% The logistic model is appealingly simple and adequate for some
% situations, but it is far too generic to capture other phenomena. Another
% way to limit growth is the *Gompertz model*, in which, for example,
GompertzPCRate = @(P) 0.5*log(6./P)/log(6/0.2);

%%
% Compared to the logistic model, the Gompertz model shuts down growth more
% rapidly until $P$ gets close to the carrying capacity.
plot( chebfun(logisticPCRate,[0.2 6]) ), hold on
plot( chebfun(GompertzPCRate,[0.2 6]) )
xlabel('P'), ylabel('P''/P')

%%
% The solutions reflect this difference, of course.
N = chebop( @(t,P) diff(P) - P.*GompertzPCRate(P), [0 25] );
N.lbc = 0.2;
GompertzResult = N\0;
clf, plot( [exponentialResult, logisticResult, GompertzResult] )
ylim([0 8])
xlabel('t'), ylabel('P(t)'), legend('exp','logistic','Gompertz')

%%
% The Gompertz model has been recognized for some time as a reasonable model
% for some tumors [1-2]. However, it too is considered inadequate for many
% realistic uses. At least some of the attention paid to the logistic and
% Gompertz models has to be attributed to the fact that they have analytic
% solutions.

%%
% References:
%
% 1. Laird, A. K. Dynamics of growih in tumors and in normal organisms.
% Nati. Cancer Inst. Monogr. 30: 15-28, 1969.
%
% 2. Winsor, C. P. The Gompertz curve as a growih curve. Proc. Natl. Acad. 
% Sci. USA 18: 1-7, 1932.

