%% DEMO OF PIECEWISE OPERATORS
% Nick hale, 9 November 2010

%%
% (Chebfun Example ode/PiecewiseLinopDemo.m)

%%
% Here we demonstrate piecewise operators (incl. boundary conditions), 
% and how the systems constructor goes about solving them.

clc, clear, close all

%%
% Define a domain, and an operator with a jump at x = 0.
[d x] = domain(-1,1);
A = -diff(d,2) + diag(sign(x))

%%
% Although d has no breakpoint, diag inherits the break from sign, which
% is then inherited by A in plus.

%%
% What does it look like when we evaluate these piecewise linops?
% Well without boundary conditions, with just have two independant blocks.
A([5 4],'nobc')

%%
% But we choose boundary conditions to enforce continuity of derivatives
% up to the differential order of the operator.
A([5 4],'bc')

%%
% So here we see that continuitity of the solution and it's derivative are
% enforced across the break at x = 0.

%%
% However, we still need to apply some boundary conditions of our own to 
% the operator. Let's choose dirichlet for simplicity.
B = A & 'dirichlet';
B([5 4],'bc')

%%
% We're now almost in a position to start solving piecewise ODEs. However,
% the standard constructor won't quite cut it here, as when it constructs
% on a domain such as [-1 0 1], the 2 subdomains are treated independantly.
% By wrapping the domain as a cell, we force the use of the systems
% constructor which doesn't suffer from this. 
%
% Note also that myfun returns a cell. This is to allow the systems
% constructor to be as general as possible, however, is it likely that it
% will ever be used for anything but this? In which case we might want to
% simplify the notation? 

myfun = @(x,N,bks) B(N{:},'bc')\[ones(sum(N{:})-4,1) ; zeros(4,1)];
u = chebfun(myfun,{[-1 0 1]},'eps',1e-9);

%%
% An alternative notation is
%   u = chebfun(myfun,[-1 0 1],'eps',1e-9,'sys',1);

%%
% Of course most users won't even see things at this level - they'll just 
% be calling backslash!
v = B\1;

%%
% which we see does much the same as above.
plot(u,'b',v,'--r','LineWidth',1.6)



