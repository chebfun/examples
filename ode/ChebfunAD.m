%% Automatic differentiation in Chebfun
% Asgeir Birkisson, 16th November 2011

%%
% (Chebfun example ode/ChebfunAD.m)
% [Tags: #automaticdifferention, #AD]

%%
% The key to many capabilities of Chebfun, such as solution of nonlinear
% differential equations and detecting whether a given operator is linear
% or not is based on Chebfun's ability to compute derivatives of operators.
% The derivative computation is achived via automatic differentiation (AD),
% which differs from both traditional numerical and symbolic
% differentiation in that it's both accurate and fast. For a beginners
% introduction to AD, see [1].

%%
% In short, AD in Chebfun works by storing information about derivatives in
% every chebfun created. This information is stored in the .jacobian fields
% of chebfuns. By using the chain rule from calculus, the information is
% then used to accumulate partial derivatives and finally returning the
% derivative of one chebfun with respect to another.

%%
% As chebfuns are built up by operations on previously defined functions,
% an evaluation tree is created, giving the user an idea how the chain rule
% will eventually be used to compute the derivate. Chebfun offers the
% option to display the AD information, both in graphical and text form.
% This short example describes how Chebfun can be used to visualise how the
% chain rule will be used in action.

%% A simple example
%
% As often in Chebfun computation, we start by creating the linear function
% x on our domain of interest
x = chebfun('x');
%%
% We then build up more complicated functions from that function
u = sin(x);
%%
% To look at the derivative information now stored in u, we create a new
% variable whose value will be the .jacobian field of u:
ujac = u.jacobian
%%
% Looking at the output, we see that ujac contains a list of instructions
% on how to compile the derivative -- some are used for linearity detection
% (not of interest to us at the moment) and other for computing the
% derivative itself.

%%
% Most importantly, we see that that one part of the derivative involves
% the cosine function, and that the other part is a call to diff. This
% demonstrates the chain rule in action, as this corresponds to the
% derivative one expects to obtain when differentation the sin function in
% ordinary calculus. Further down the output, we see a mention to an "empty
% anon", this indicates the bottom of the evaluation tree, i.e. the
% independent variable x.

%% A slightly more complicated example
%
% We know introduce a new chebfun in the mix, to see how the derivatives
% look as we build bigger evaluation trees
v = cos(u) + x.*u;
%% 
% The derivative information of v takes up more lines than in the previous
% example
vjac = v.jacobian
%%
% It is easy to imagine the tree getting very large for more involved
% computation. Hence, we introduce the plot method for AD information which
% plots the tree in a graphical way:
plot(vjac)
%%
% The tree in the figure shows how the final chebfun v is composed by
% performing various operations on the chebfun we started with, x. Note
% that by clicking on the nodes of the tree, the AD information shown
% before is displayed.

%% A nice evaluation tree
%
% Finally, we create the complicated chebfun f via
f = exp(x.^2) + cos(x).*log(2+x) + diff(tanh(x.*u));
%%
% Here, printing the AD information to the console would be difficult to
% grasp, but we get a nice tree if we plot it graphically:
plot(f.jacobian)

%%
% References:
% [1] http://en.wikipedia.org/wiki/Automatic_differentiation.