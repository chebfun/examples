%% Best polynomial approximation in the L1 norm 
% Yuji Nakatsukasa and Alex Townsend, July 2019

%%
% (Chebfun example approx/polyfitL1.m)
% [Tags: #polynomial, #best, #Newton]

%% Best polynomial approximation in the $L_\infty$ norm 
% Given a continuous real-valued function $f$ on $[a,b]$, finding the best
% polynomial approximant to $f$ in the $L_\infty$-norm 
% $\mbox{min}_{p\in\mathcal{P}_n}\|f-p\|_{L_{\infty}}$
% is known as the minimax (sometimes Chebyshev
% approximation) problem. It is computed in Chebfun by the command
% `minimax()`, and the error $f-p_\infty$ exhibits the beautiful
% equioscillation phenomena. Let's revisit the Chebfun example
% {\tt https://www.chebfun.org/examples/approx/ResolutionWiggly.html} and 
% compute its best polynomial approximant. 

dom = [0 14]; format compact
f = chebfun(@(x) sin(x).^2 + sin(x.^2), dom);
LW = 'linewidth'; lw = 1.2;
plot(f, LW, lw, LW, lw), ylim([-2.5 2.5])
hold on

deg = 80; 
pinf = minimax(f, deg, 'maxiter', 100); 
plot(pinf, LW, lw)

%% Best polynomial approximation in the $L_2$ norm 
% The best polynomial approximant to $f$ in the $L_2$-norm is algorithmically
% easier to compute as it is the orthogonal projection of $f$ onto the space
% of polynomials of degree $\leq n$. In Chebfun, the `polyfit()` command 
% computes this best polynomial approximant: 

p2 = polyfit(f, deg); 
hold on 
plot(p2, LW, lw)

%% Best polynomial approximation in the $L_1$ norm 
% Recently, Chebfun added a command `polyfitL1()` to compute the best polynomial
% approximant to $f$ in the $L_1$-norm. (See Pinkus' book for a survey of
% results about best $L_1$ polynomial approximants~[2].) Compressing sensing 
% has made the $L_1$-norm an important tool in signal processing as it can 
% promote sparsity in the solution or residue. A Newton-based algorithm proposed by 
% Watson~[4] is known to converge to the best polynomial approximant in the 
% $L_1$-norm, under some assumptions. 

p1 = polyfitL1(f, deg);
plot(p1, LW, lw)
legend('f', 'L_\infty', 'L_2', 'L_1', 'Location', 'Best')
set(gca, 'FontSize', 14)

%%
% To highlight the difference between the best $L_\infty$, $L_2$, and $L_1$
% polynomial approximants of $f$, we plot the errors together.

clf
plot(pinf-f, LW, lw), hold on
plot(p2-f, LW, lw)
plot(p1-f, LW, lw) 
title('Error plot')
legend('L_\infty', 'L_2', 'L_1', 'Location', 'NorthWest')
set(gca, 'FontSize', 14);

%%
% Are these plots interesting? In particular, do they differ significantly?
% To gain some insight we zoom in to the left-part of the figure:
xlim([0,7])

%%
% One can see that the error with $L_1$ is much smaller in this region. We 
% say that the $L_1$ error is more localized. 

%% Another example
% Let's do one more example, motivated by [3, Appendix, Myth 3]. Let's try to 
% approximate $|x-1/4|$ on $[-1,1]$ with a polynomial of degree
% $\leq 80$. Below, we plot the error functions with best $L_\infty$,
% $L_2$, and $L_1$ polynomial approximants: 

f = chebfun(@(x) abs(x-1/4),'splitting','on');

deg = 80; 
pinf = minimax(f, deg, 'maxiter', 100);
p2 = polyfit(f, deg); 
p1 = polyfitL1(f, deg);

clf 
plot( f-pinf, LW, lw ), hold on
plot( f-p2, LW, lw)
plot( f-p1, LW, lw)
title('Error plot')
legend('L_\infty', 'L_2', 'L_1', 'Location', 'SouthWest');
set(gca, 'FontSize', 14);

%%
% Let's zoom in to the left half portion of the figure. 
xlim([-1,-0.4])

%%
% Again, we see that the best $L_1$ polynomial approximant has a far more 
% localized error. This is actually a typical phenomenon that is explained 
% in~[1]. 

%% A short word on the underlying algorithm for polyfitL1()
% In [1], it is advocated that Watson's algorithm should be used in 
% conjunction with linear programming problems and a refinement step. These
% additional algorithmic details can significantly speed up Watson's algorithm. 
% Since the linear programming commands in MATLAB are in a commerical 
% toolbox, we have decided to avoid them so that all Chebfun users can 
% execute the `polyfitL1()` command. 

%% References 
% [1] Y. Nakatsukasa and A. Townsend, Error localization of best L1
% polynomial approximants, arXiv:1902.02664. 
% 
% [2] A. M. Pinkus, On L1-approximation. Vol. 93. Cambridge University 
% Press, 1989.
%
% [3] L. N. Trefethen, Approximation Theory and Approximation Practice,
% SIAM, 2012. 
% 
% [4] G. A. Watson. An algorithm for linear L1 approximation of continuous functions. IMA J.
% Numer. Anal., 1(2):157-167, 1981.