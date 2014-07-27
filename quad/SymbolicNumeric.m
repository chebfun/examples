%% Symbolic and numerical integration
% Nick Trefethen, July 2014

%%
% (Chebfun example quad/SymbolicNumeric.m)
% [Tags: #quadrature, #Mathematica, #SUM]

%%
% Mathematica can do extraordinary things with symbolic
% integration.  One way to see it in action is to go
% to the Wolfram Mathematica Online Integrator and
% click on "Random Example".

%%
% For example, consider the function
%
% $$ f(x) = \log(2+x)^3 \log(3+x) x^3 $$
%
% The Online Integrator quickly delivers the following exact indefinite
% integral:

%%
% $$ (-759744 - 558290 x + 17705 x  - 1050 x  + 54 x  + \\\\
%   910528 \text{Log}[2 + x] + 400008 x \text{Log}[2 + x] - \\\\
%   22836 x^2  \text{Log}[2 + x] + 2072 x^3  \text{Log}[2 + x] - \\\\
%   162 x^4  \text{Log}[2 + x] - 302016 \text{Log}[2 + x]^2  - \\\\
%   118800 x \text{Log}[2 + x]^2  + 11880 x^2  \text{Log}[2 + x]^2  - \\\\
%   1680 x^3  \text{Log}[2 + x]^2  + 216 x^4  \text{Log}[2 + x]^2  + \\\\
%   48384 \text{Log}[2 + x]^3  + 15552 x \text{Log}[2 + x]^3  - \\\\
%   2592 x^2  \text{Log}[2 + x]^3  + 576 x^3  \text{Log}[2 + x]^3  - \\\\
%   144 x^4  \text{Log}[2 + x]^3  + 309078 \text{Log}[3 + x] + \\\\
%   79680 x \text{Log}[3 + x] - 5520 x^2  \text{Log}[3 + x] + \\\\
%   592 x^3  \text{Log}[3 + x] - 54 x^4  \text{Log}[3 + x] - \\\\
%   293976 \text{Log}[2 + x] \text{Log}[3 + x] - \\\\
%   57600 x \text{Log}[2 + x] \text{Log}[3 + x] + \\\\
%   7488 x^2  \text{Log}[2 + x] \text{Log}[3 + x] - \\\\
%   1344 x^3  \text{Log}[2 + x] \text{Log}[3 + x] + \\\\
%   216 x^4  \text{Log}[2 + x] \text{Log}[3 + x] + \\\\
%   138672 \text{Log}[2 + x]^2  \text{Log}[3 + x] + \\\\
%   13824 x \text{Log}[2 + x]^2  \text{Log}[3 + x] - \\\\
%   3456 x^2  \text{Log}[2 + x]^2  \text{Log}[3 + x] + \\\\
%   1152 x^3  \text{Log}[2 + x]^2  \text{Log}[3 + x] - \\\\
%   432 x^4  \text{Log}[2 + x]^2  \text{Log}[3 + x] - \\\\
%   46656 \text{Log}[2 + x]^3  \text{Log}[3 + x] + \\\\
%   576 x^4  \text{Log}[2 + x]^3  \text{Log}[3 + x] - \\\\
%   24 (5609 - 6756 \text{Log}[2 + x] + 4680 \text{Log}[2 + x]^2 ) \\\\
%    \text{PolyLog}[2, -2 - x] + \\\\
%   288 (-563 + 780 \text{Log}[2 + x]) \text{PolyLog}[3, -2 - x] - \\\\
%   224640 \text{PolyLog}[4, -2 - x]) / 2304 $$

%%
% In Chebfun, more prosaically, we could do this:
LW = 'LineWidth'; CO = 'Color'; FS = 'FontSize';
f = chebfun(@(x) log(2+x).^3.*log(3+x).*x.^3);
fi = cumsum(f)
plot(fi,LW,2.2,CO,[0 .7 0]), grid on
title('symbolically integrable',FS,14)

%%
% These two results are so utterly different! -- and each would
% be superior in some applications.

%%
% The definite integral from $-1$ to $1$ could be computed
% in Chebfun like this,
sum(f)

%%
% or like this,
fi(1) - fi(-1)

%%
% If I ask WolframAlpha for the definite integral,
% it gives six numerical digits, $0.364264$, so perhaps it
% is bypassing the symbolic solution.

%%
% An example like this highlights the combinatorial complexity
% that can arise in symbolic computing.  Of course, sometimes
% no symbolic answer is available at all.
% If $f$ is changed to
%
% $$ g(x) = \log(2+x)^3 \log(3+x)^2 x^3, $$
%
% then the Online Calculator reports: "Mathematica could not
% find a formula for your integral.  Most likely this means
% that no formula exists."  For Chebfun, on the other
% hand, it makes no difference:
g = chebfun(@(x) log(2+x).^3.*log(3+x).^2.*x.^3);
gi = cumsum(g)
plot(gi,LW,2.2,CO,[.7 0 .7]), grid on
title('not symbolically integrable',FS,14)

%% References
%
% 1.  L. N. Trefethen, Computing numerically with functions
%     instead of numbers, _Communications of the ACM_, to appear.
