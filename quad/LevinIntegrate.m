%% An oscillatory integral
% Sheehan Olver, December 2010

%%
% (Chebfun example quad/LevinIntegrate.m)
% [Tags: #quadrature, #ODE]

%
% This example computes the highly oscillatory integral of
%
%       f * exp( 1i * w * g ),
%
% over (0,1) using the Levin method [1]. This method computes the integral
% by rewriting it as an ODE
%
%       u' + 1i * w * g' u = f,
%
% so that the indefinite integral of f * exp( 1i * w * g ) is 
%
%       u * exp( 1i * w * g ).


%%
% We use as an example
%
%       f = 1 / ( x + 2 );
%       g = cos( x - 2 );
%       w = 100000;

d = domain(0,1);
x = chebfun(@(x) x, d);
f = 1./(x+2);
g = cos(x-2);
D = diff(d);

%%
% Here is are plots of this integrand, with w = 100, in complex space
w = 100;
LW = 'LineWidth'; lw = 1.6; FS = 'FontSize'; fs = 14;
plot(f.*exp(1i*w*g), LW, lw), axis equal
title('Complex plot of integrand',FS,fs)
%%
% and of just the real part
plot(real(f.*exp(1i*w*g)), LW, lw), axis equal
title('Real part of integrand',FS,fs)

%%
% The Levin method will be accurate for large and small w, and the time
% taken is independent of w. Here we take a reasonably large value of w.
w = 100000;

%%
% Start timing
tic

%%
% Construct the operator L
L = D + 1i*w*diag(diff(g));

%%
% Since we have no boundary conditions, we fake that the derivative order
% is zero.
L = set(L, 'difforder', 0);

%%
% From asymptotic analysis, we know that there exists a solution to the
% equation which is non-oscillatory, though we do not know what initial
% condition it satisfies.  Thus we find a particular solution to this
% equation with no boundary conditions.

u = L \ f;

%%
% Because L is a differential operator with derivative order 1, \ expects
% it to be given a boundary condition, which is why the warning message is
% displayed. However, this doesn't cause any problems: though there are,
% in fact, a family of solutions to the ODE without boundary conditions
% due to the kernel
%
%     exp(- 1i * w * g),
%
% it does not actually matter which particular solution is computed.
% Non-uniqueness is also not an issue: \ in matlab is least squares, hence
% does not require uniqueness. The existence of a non-oscillatory solution
% ensures that \ converges to a u with length independent of w.
%
% One could prevent the warning by applying a boundary condition consistent
% with the rest of the system, that is 
%  L.lbc = {L(1,:),f(0)};

%%
% Now we evaluate the antiderivative at the endpoints to obtain the
% integral.

u(1).*exp(1i.*w.*g(1)) - u(0).*exp(1i.*w.*g(0))

toc

  
%%
% Here is a way to compute the integral using Clenshaw--Curtis quadrature.
% As w becomes large, this takes an increasingly long time as the
% oscillations must be resolved.

tic
sum( f.*exp(1i.*w.*g) )
toc

%%
% References:
%
% [1] Levin, D., Procedures for computing one and two-dimensional integrals
% of functions with rapid irregular oscillations, Maths Comp.,  38 (1982) 531--538
