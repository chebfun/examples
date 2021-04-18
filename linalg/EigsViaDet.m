%% Computing eigenvalues by sampling the determinant
% Jared Aurentz and Nick Trefethen, November 2014

%%
% (Chebfun Example linalg/EigsViaDet.m)
% [Tags: #linearalgebra, #eigenvalues, #determinant]

function EigsViaDet()
%%
% The eigenvalues of a matrix $A$ are the roots of the 
% determinant function, $f(x) = \det(xI-A)$.  If $A$ is
% real symmetric and tridiagonal and of dimension $N$,
% then $f(x)$ can be computed in $O(N)$ operations by
% the method known as Sturm sequences, described in
% many texts such as on p. 229 of [2] or p. 423 of [3].
% Let $d_k$ denote the determinant of the upper-left $k\times k$
% block of $A$, and let the diagonal and superdiagonal
% entries of $A$ be $a_k$ and $b_k$, respectively.
% Then the key observation, easily derived, is
% that the determinants satisfy a 3-term recurrence
% relation,
%
% $$ d_{n+1} = a_{n+1} d_n - b_n^2 d_{n-1} . $$
%
% It is salutary to note that we can easily 
% vectorize this recurrence, which means that 
% Chebfun can use it very efficiently to construct
% chebfuns corresponding to $f(x)$ over a prescribed
% interval.

%%
% Here is our function for evaluating $f(x)$, which
% we will call `fdet`.
     function fdet = fdet(x,a,b,N)
         dold = ones(size(x));
         d = x-a(1);
         for k = 1:N-1
             dnew = (x-a(k+1)).*d - b(k)^2.*dold;
             dold = d; d = dnew;
         end
         fdet = d;
     end
            
%%
% OK, let's try it.  Here is a matrix whose eigenvalues
% lie roughly in the interval $[-5,5]$:
tic
N = 100;
rng(2)
a = 10*rand(N,1)-5;
b = randn(N-1,1);
A = spdiags([[b;0] a [0;b]],-1:1,N,N);

%%
% Here, computed the usual way, are
% the "exact" eigenvalues in the interval $[-1,1]$:
format long
e = eig(full(A)); e_exact = sort(e(abs(e)<=1))

%%
% Here we make a chebfun of the determinant function:
c = chebfun(@(x) fdet(x,a,b,N),[-1,1]);
plot(c), grid on
xlabel('x')
title('det(xI-A) as a chebfun')

%%
% Now we compute its roots and compare them
% with the true eigenvalues.
e_inexact = roots(c);
disp('         exact              inexact            difference')
disp([e_exact e_inexact e_exact-e_inexact])

%%
% Is this good agreement?  Well things
% look pretty good, but for many of the eigenvalues we are losing up to 
% five digits of accuracy, and in fact, this method faces difficulties and
% would quickly fail for larger values of $N$.
% A plot of the absolute value of |c| on a log scale
% gives an indication of what is going on.
semilogy(abs(c)), ylim([1e22 1e32]), grid on
xlabel('x')
title('|det(xI-A)| on a log scale')

%%
% The first thing we note in this figure is that
% the scale of the data is a long way from $1$.  This has
% something to do with the scaling of the problem to the
% interval $[-5,5]$, and could be alleviated to some
% extent by a rescaling.  It could only be alleviated
% partially, however, for the more fundamental problem is the
% exponential variation of scales across the interval, a phenomenon
% associated with the subject of potential theory [1].
% This is a mathematical fact about the determinant function.
% Wilkinson pointed out that in fact the determinant function
% can be computed with high relative accuracy, despite the bad
% scaling [3, p. 228], so the problem in our method is not its reliance
% on $det(xI-A)$.  Rather, it is in making a chebfun representation
% of this function over a broad interval.

%%
% To confirm this, note how much better the accuracy becomes if
% we restrict attention to $[-1,0]$:
e_exact = sort(e(e<0 & abs(e)<1));
c = chebfun(@(x) fdet(x,a,b,N),[-1,0]);
plot(c), grid on, ylim([-5e26 1e27])
xlabel('x')
title('det(xI-A) on a smaller interval')
e_inexact = roots(c);
disp('         exact              inexact            difference')
size(e_exact), size(e_inexact)
disp([e_exact e_inexact e_exact-e_inexact])

%%
% Another amusing approach is to use Chebfun's edge detector
% to count eigenvalues!  The accuracy is magnificent, showing
% that Chebfun's edge detector is not thrown off by bad scaling.
c2 = chebfun(@(x) sign(fdet(x,a,b,N)),[-1,1],'splitting','on','minSamples',100);
plot(c2,'jumpline','-'), grid on, ylim([-1.4 1.4]);
e_edgedetect = roots(c2);
hold on, plot(e_edgedetect,0*e_edgedetect,'.r'), hold off
disp('         exact        via edge detection      difference')
e_exact = sort(e(abs(e)<=1));
size(e_exact), size(e_edgedetect)
disp([e_exact e_edgedetect e_exact-e_edgedetect])

%%
% Here is the total time for this Example:
toc
end

%% References
%
% 1. L. N. Trefethen, _Approximation Theory and Approximation Practice_, SIAM,
%    2013.
%
% 2. L. N. Trefethen and D. Bau, III, _Numerical Linear Algebra_, SIAM, 1997.
%
% 3. J. H. Wilkinson, _The Algebraic Eigenvalue Problem_, Clarendon Press,
%    1965.
