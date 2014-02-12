%% Jump Conditions in BVPS
% Nick Hale, 25th November 2011

%%
% (Chebfun Example ode/JumpConditions.m)
% [Tags: #linearODE, #jumpconditions]

%%
% Chebfun has recently added support for jump conditions in solutions to
% differential equations. Here we'll demonstrate this on a few toy examples.

%%
% Let's start with a simple linear ODE

N = chebop(@(x,u) 1e-2*diff(u,2) + sin(x).*u);

%%
% and add some dirichlet boundary conditions

N.lbc = 1; N.rbc = 1;

%%
% Now let's solve this as it is (i.e., without jumps) for reference.

u1 = N\0
plot(u1)

%%
% Suppose we want to add a jump condition at the origin. To do this we'd
% use the .bc field as follows

N.bc = @(x,u) feval(u,0,'right') - feval(u,0,'left') - 1;
u2 = N\0
plot(u2)

%%
% The above notation is a bit clunky, and the syntax 'jump' can be used
% instead. For example, we get the same result as above with

N.bc = @(x,u) jump(u,0) - 1;
u3 = N\0
norm(u3-u2)

%%
% We can quickly around and make a jump appear instead in the derivative

N.bc = @(x,u) jump(diff(u),0) + 2*pi;
u4 = N\0
plot(u4)

%%
% Or go crazy and introduce multiple jumps!

N.bc = @(x,u) [jump(u,-.8:.2:.8) - (-.8:.2:.8),...
               jump(diff(u),-.8:.2:.8) + (-.8:.2:.8)];
u5 = N\0
plot(u5)

%%
% I'd be very interested to hear if you have any practical problems which
% require these kinds of conditions!
