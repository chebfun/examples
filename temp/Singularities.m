%% Types of singularities
% Hrothgar, 12th November 2013

%%
% (Chebfun example complex/Singularities.m)
% [Tags: #complex, #phase portraits, #Chebfun2]

function Singularities
LW = 'linewidth'; FS = 'fontsize'; MS = 'markersize';

%%
% A singularity $z_0$ of a function $f$ in the complex plane
% is classified as removable, pole of order $n$, or essential
% depending on the coefficients of the main part of the
% function's Laurent series expansion at $z_0$,
% $$ f(z) = \sum_{k=-\infty}^\infty c_k (z - z_0)^k. $$
% Specificially, we say $z_0$ is [1]
%
% 1. a removable singularity if $c_k = 0$ for all $k < 0$,
%
% 2. a pole of order $n$ if $c_{-n} \neq 0$ and $c_k = 0$ for all ...
% $k < -n < 0$,
%
% 3. an essential singulatiry if $c_k \neq 0$ for infinitely many ...
% negative $k$.
%
% Equivalently, we can characterize a singularity $z_0$ of $f$ by
% the limits of $f$ and $\frac1f$ at $z_0$:
%
% 1. If both $\lim_{z\rightarrow z_0} f(z)$ and $\lim_{z\rightarrow z_0} \frac1{f(z)}$
% exist, then $z_0$ is a removable singularity of $f$.
%
% 2. If $\lim_{z\rightarrow z_0} f(z)$ does not exist but
% $\lim_{z\rightarrow z_0} \frac1{f(z)}$ does exist, then $z_0$ is a pole of $f$.
%
% 3. If neither $\lim_{z\rightarrow z_0} f(z)$ nor $\lim_{z\rightarrow z_0} \frac1{f(z)}$
% exists, then $z_0$ is an essential singularity of $f$.

%% A removable singularity
% Consider the function $f(z) = \mathrm{sinc}(z) = \frac{\sin(z)}{z}$. Strictly
% speaking, $f$ is not defined at $z = 0$. However, expanding $f$ as
% a Laurent series reveals that there are no nonzero coefficients
% for negative powers of $z$:
% $$ f(z) = \frac1z \sum_{k=0}^\infty \frac{(-1)^k z^{2k+1}}{(2k+1)!} $$
% $$ \phantom(f(z)) = \sum_{k=0}^\infty \frac{(-1)^k z^{2k}}{(2k+1)!} $$
% $$ \phantom(f(z)) = 1 - \frac{z^2}{3!} + \frac{z^4}{4!} - \frac{z^5}{6!} + \cdots. $$
% The singularity at $z = 0$ can be removed simply by defining
% $f(0) = c_0 = 1$.

f = @(z) sinc(z/pi); % sinc in MATLAB is defined as sin(pi z)/(pi z)
removable = chebfun2(f, 1.5*pi*[-1 1 -1 1]);
plot(removable)

% Indeed, the phase portrait of $f$ looks clean around
% the origin, and there is no singularity after all. The points
% $\pm pi$ stand out in the plot, but the are zeros, not singularities,
% of $f$. In a phase portrait the difference is that the colors of
% zeros and poles wind in opposite directions.

%% Poles
% Contrary to removable singularities, we can read off
% singularitities of the form $\frac1{z^n}$ from a phase portrait.
% In order for Chebfun2 to handle the poles, we will "smash"
% the function á la Nick Trefethen's earlier example "Phase portraits
% for functions with poles" -- that is, we will plot a smooth
% function with the same phase as the one we're interested in.

function g = smash(f)
    g = f./(1 + abs(f).^2); % smooth function with same phase as f
    g(isnan(g)) = 0;        % give 0 rather than NaN at poles
end

% Now we create a function with poles of different orders at
% the points ${\pm 1, \pm \mathrm{i}}$.

g = @(z) (z-1).^-1 .* (z-1i).^-2 .* (z+1).^-3 .* (z+1i).^-4;
poles = chebfun2(@(z) smash(g(z)), 2*[-1 1 -1 1]);
plot(poles)

%%
% The order of each pole is equal to the number of times each
% color appears when winding once around the pole. At $z = 1$, each
% color appears once, indicating that $z = 1$ is a simple pole of $f$.
% On the other hand, at $z = -1$, each color appears three times,
% indicating that $z = -1$ is a pole of order 3.

%% An essential singularity
% The function $h(z) = \mathrm{e}^\frac1z$ has an essential singularity at
% the origin. As mentioned above, this means its Laurent expansion
% at $z=0$ has infinitely many negative terms. In particular,
% $$ h(z) = \sum_{k=-\infty}^0 \frac{(-z)^k}{(-k)!}.$$
% This unwieldy singularity cannot be captured in full by Chebfun2,
% but we can peek at it from the side by squashing the complex
% plane with a transformation $z \mapsto z^\frac9{10}$.

h1 = @(z) exp(-1./(z.^.9));
essential1 = chebfun2(h1, .5*[0 1 -.5 .5]);
plot(essential1)

% The essential singularity can be thought of as a pole of order infinity.
% Winding around at an infinitesimal distance from the origin, each color
% appears infinitely many times.
% 
% [statement of Big Picard Theorem]

%%
% Finally, consider the function $h(z) = \sin(\frac1z)$, which has a limit
% point of singularities at the origin. [...]

h2 = @(z) sin(1./(z.^.8));
essential2 = chebfun2(@(z) smash(h2(z)), .5*[0.02 1 -.5 .5]);
plot(essential2)

%%
% References:
%
% [1] Wegert, Elias. Visual complex functions. Vol. 1. Springer, 2012.
%
% [2] M. J. Ablowitz and A. S. Fokas. Complex variables: introduction 
%     and applications. Cambridge University Press, 1997.
%
end
