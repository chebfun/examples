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
%
% (-759744 - 558290 x + 17705 x  - 1050 x  + 54 x  + 
%
%   910528 Log[2 + x] + 400008 x Log[2 + x] - 
%
%          2                    3
%   22836 x  Log[2 + x] + 2072 x  Log[2 + x] - 
%
%        4                               2
%   162 x  Log[2 + x] - 302016 Log[2 + x]  - 
%
%                      2          2           2
%   118800 x Log[2 + x]  + 11880 x  Log[2 + x]  - 
%
%         3           2        4           2
%   1680 x  Log[2 + x]  + 216 x  Log[2 + x]  + 
%
%                   3                     3
%   48384 Log[2 + x]  + 15552 x Log[2 + x]  - 
%
%         2           3        3           3
%   2592 x  Log[2 + x]  + 576 x  Log[2 + x]  - 
%
%        4           3
%   144 x  Log[2 + x]  + 309078 Log[3 + x] + 
%
%                              2
%   79680 x Log[3 + x] - 5520 x  Log[3 + x] + 
%
%        3                  4
%   592 x  Log[3 + x] - 54 x  Log[3 + x] - 
%
%   293976 Log[2 + x] Log[3 + x] - 
%
%   57600 x Log[2 + x] Log[3 + x] + 
%
%         2
%   7488 x  Log[2 + x] Log[3 + x] - 
%
%         3
%   1344 x  Log[2 + x] Log[3 + x] + 
%
%        4
%   216 x  Log[2 + x] Log[3 + x] + 
%
%                    2
%   138672 Log[2 + x]  Log[3 + x] + 
%
%                     2
%   13824 x Log[2 + x]  Log[3 + x] - 
%
%         2           2
%   3456 x  Log[2 + x]  Log[3 + x] + 
%
%         3           2
%   1152 x  Log[2 + x]  Log[3 + x] - 
%
%        4           2
%   432 x  Log[2 + x]  Log[3 + x] - 
%
%                   3
%   46656 Log[2 + x]  Log[3 + x] + 
%
%        4           3
%   576 x  Log[2 + x]  Log[3 + x] - 
%
%                                               2
%   24 (5609 - 6756 Log[2 + x] + 4680 Log[2 + x] ) 
%
%    PolyLog[2, -2 - x] + 
%
%   288 (-563 + 780 Log[2 + x]) PolyLog[3, -2 - x] - 
%
%   224640 PolyLog[4, -2 - x]) / 2304

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

fi = cumsum(f)

%%
% Reference:
% 
% 1.  L. N. Trefethen, Computing numerically with functions
%     instead of numbers, _Communications of the ACM_, to appear.
