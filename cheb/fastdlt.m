%% The fast discrete Legendre transform 
% Nick Hale and Alex Townsend, April 2015 

%% 
% (Chebfun example cheb/fastdlt.m)
% [Tags: #Legendre, #transform, #nonuniform]

FS = 'fontsize'; 
LW = 'linewidth'; 

%% The forward transform 
% The discrete Legendre transform (DLT), or finite Legendre transform, evaluates
% a Legendre series expansion at Legendre nodes on $[-1,1]$, i.e., 
% 
%  $$ f(x_k^{leg}) = \sum_{n=0}^{N-1} c_n^{leg} P_n( x_k^{leg} ), \qquad 0\leq k\leq N-1 $$
% 
% It is considered an awkward task because Legendre nodes are 
% non-uniform and Legendre polynomials have no explicit closed-form 
% expression. Therefore, FFT-based algorithms do not
% immediately apply and until recently the transform required $O(N^2)$ 
% operations. (See [3] for one of the earliest fast DLT algorithms.)

%% 
% In [2] we describe an algorithm to compute the transform in 
% $O(N(log N)^2/loglogN)$ operations, which is implemented in the dlt 
% command in Chebfun.  This allows us to compute 
% the transform when N is 10000 or 100000, or a million. Here it is in 
% action: 

cleg = randn( 1e4, 1); 
tic, fleg = chebfun.dlt( cleg ); toc

%% 
% The transform computed above is still based on the FFT, but with a
% handful of tricks and approximations to make it applicable to the
% DLT. Two main facts are exploited.

%% A Legendre polynomial can be related to a cosine
% An asymptotic expansion of a Legendre polynomial of high degree shows
% that $P_n$ is not too far away from a cosine (after a change of
% variables and scaling). The precise statement is, as
% $n\rightarrow\infty$, 
% 
%   $$ \sqrt{\sin(\theta)}P_n(\cos(\theta)) \sim  \cos( (n+1/2)\theta +
%   (n-1/4)\pi ). $$
%
% A signal processing engineer would verify this as follows: 

P = chebfun( leg2cheb([zeros(1e4-1,1) ; 1]), 'coeffs' ); % Legendre polynomial
theta = linspace(0,2*pi,4e4);                      % time samples
modifiedSignal = sqrt(sin(theta)).*P(cos(theta));  % modify signal 
modes = abs( fft( modifiedSignal ) );              % amplitude 
plot( modes(1:2e4), LW, 2 ), hold on,              % freq-domain plot
xlabel('Frequency bins',FS,14) 
ylabel('Magitude',FS,14) 
title('Frequency analysis of modified Legendre polynomial',FS,10)
set(gca,FS,14)

%% 
% The engineer would go on to say that since only a handful of modes near
% 1e4 are excited, the signal $\sqrt{\sin(\theta)}P_n(\cos(\theta))$ can be
% well-approximated by a handful of sinusoidal waves. 

%%
% This is one fact that is exploited in our fast algorithm for computing 
% a discrete Legendre transform. 

%% Legendre nodes are nearly uniform on the unit circle
% Another fact is that Legendre nodes are well-approximated by Chebyshev
% nodes (Chebyshev points of the first kind), as we can see below:

NN = floor(logspace(1,4,50));
j = 1; 
for N = NN 
    xleg = acos(legpts( N )); 
    xcheb = acos(chebpts( N, 1));
    maxdiff(j) = norm( xleg - xcheb, inf ); j = j  + 1; 
end

loglog( NN, maxdiff, LW, 2 ), hold on 
loglog( NN, 0.83845./NN, 'k--', LW, 2)
legend('||x^{leg} - x^{cheb}||_\infty', 'Theoretical bound') 
xlabel('N', FS, 14), ylabel('Max abs diff',FS, 14)

%% 
% These two facts mean that the DLT can be carefully
% related to a handful of discrete cosine transforms, allowing it be
% calculated via the FFT, requiring $O(N(log N)^2/loglogN)$ operations. The
% odd-looking complexity comes about because of a balancing of computation 
% costs, see [1].

%% The inverse transform 
% The inverse discrete Legendre transform takes samples of a function 
% from a Legendre nodes and computes the associated Legendre expansion
% coefficients.  

f = chebfun( @(x) 1./(1+10000*x.^2) ); 
cleg = cheb2leg( f.coeffs ); 
tic, fleg = chebfun.dlt( cleg ); toc
backToTheCoeffs = chebfun.idlt( fleg );
norm( backToTheCoeffs - cleg, inf )

%% 
% The IDLT can be related to a (transposed) DLT by a discrete orthogonality 
% relation. It has the same $O(N(log N)^2/loglogN)$ complexity. 

%% 
% The DLT and IDLT complete the Chebyshev--Legendre cycle: (In the diagram 
% below we give the Chebfun commands that compute each particular transform.) 
%
%                   -->-->--coeffs2vals-->-->--                      
%       CHEBCOEFFS  --<--<--vals2coeffs--<--<--   CHEBVALUES
%         |  |
%         |  |      
%         ^  \/       
%        cheb2leg         
%        leg2cheb       
%         |  |
%         |  |
%       LEGCOEFFS       -->-->--dlt-->-->--      LEGVALUES
%                       --<--<--idlt--<--<--
% 
%  One can now move freely between Chebyshev and Legendre modes and values with
%  fast algorithms in Chebfun. Hooary! 

%% References
% [1] N. Hale and A. Townsend, A fast, simple, and stable
% Chebyshev--Legendre transform using an asymptotic formula, SISC, 36
% (2014), pp. A148--A167. 
% [2] N. Hale and A. Townsend, A fast FFT-based discrete Legendre
% transform, in preparation. 
% [3] D. Potts, Fast algorithms for discrete polynomial transforms on
% arbitrary grids, Linear Algebra and its Applications, 336 (2003), 
% pp. 353--370. 