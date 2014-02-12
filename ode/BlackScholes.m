%% Time independent Black-Scholes with jumps
% Alex Townsend, 28th October 2011

%%
% Chebfun example ode/BlackScholes.m
% [Tags: #linearODE, #piecewise, #jumpconditions, #blackscholes]

%% 
% The Black-Scholes equation is a partial differential equation for
% modelling the price of an European option in terms of underlying equity
% prices [1]. In this example we consider the time independent
% Black-Scholes equation which is a one-dimensional ODE.

%% Good investment?
% Let's suppose you buy an European option for £50 that depends on the
% share value of Apple Inc. and the risk-free interest rate is 3%.  At the
% time you decide to sell the option an incremental tax applies so that you
% pay 20% of the price of the share rounded down to the nearest multiple of
% 10. If the underlying share is worth £1, you lose all your investment
% and when its worth £50 you will be able to sell your option for £150.

r = 1.03;  % Risk-free interest rate
vol = 1;   % Volality
tax = 0.2; % Rate of tax
taxpts = 10:10:40;
N = chebop(@(s,V) .5*vol*s.^2.*diff(V,2) + r*s.*diff(V) - r*V,[1,50]);
N.lbc = @(V) V+50;
N.rbc = @(V) V-150;
N.bc = @(V) jump(V,taxpts)+tax*feval(V,taxpts);
y=N\0;
plot(y), hold on; 
title('Profit/loss versus underlying share price','FontSize',16);
xlabel('Share Price in pounds'); ylabel('Profit');

%% Break-even point and double your money
% As a shrewd investor, you would like to know the underlying share price
% when you break-even and when you double your money.  This can be computed
% by the roots command.

fprintf('Break-even point = £%1.2f\n',roots(y));
fprintf('Double your money = £%1.2f\n',roots(y-100));

%%
% Note that you do not double your money when the underlying share price is
% £30 this is just where the solution jumps across the £100 profit mark.

%%
% Don't lose all your money! 

%% References
% [1] http://en.wikipedia.org/wiki/Black-Scholes