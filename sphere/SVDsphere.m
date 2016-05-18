%% The singular value decomposition on the sphere
% Alex Townsend and Grady Wright, May 2016 

%% 
% (Chebfun Example spherefun)
% [Tags: #sphere, #svd, #lowrank]

%% The singular value decomposition
% <latex>
% Let $S$ denote the surface of the unit sphere and 
% $f:S\rightarrow \mathbb{R}$ be a real-valued function. If $f$ 
% is $L^2$-integrable, i.e., 
% 
% $$ \| f \|_{L^2(S)}^2 = \int_{S} |f(x,y,z)|^2 dxdydz < \infty $$,
% 
% then $f$ has a singular value decomposition (SVD) in intrinsic variables 
% on the sphere. That is, 
% 
%  \begin{equation}
%  f(\lambda,\theta) = \sum_{k=1}^\infty \sigma_k u_k(\theta)v_k(\lambda),$$
%  \label{eq:SVD}
%  \end{equation}
% 
% where $\sigma_1\geq \sigma_2\geq \cdots \geq 0$ are the singular values
% of $f$ and $u_k$ and $v_k$ are the singular functions. The equality in 
% (\ref{eq:SVD}) should be understood as holding in the mean-square sense. 

% Here, the set $\{u_1,u_2,\ldots\}$ is a complete orthonormal set of 
% integrable functions with respect to the inner-product
% 
%  $$ <p, q> = \int_{0}^\pi p(\theta)q(\theta)\sin(\theta)d\theta $$,
%
% and $\{v_1,v_2,\ldots\}$ is a complete orthonormal set of integrable 
% functions with respect to standard $L^2$ inner-product on $[-\pi,\pi]$. 
% 
% </latex>

%% 
% In Spherefun the |svd| command computes the singular value decomposition
% of a spherefun object. Below, we take a function and computes its SVD. 
% To check the accuracy of the SVD we are required to first convert to a 
% chebfun2 object.

f = spherefun( @(x,y,z) cos(x.*y.*z) ); 
[U, S, V] = svd( f ); 
norm( chebfun2(f,[-pi,pi,0,pi]) - U*S*V')

%% 
% Mathematically, most functions on the sphere are of infinite rank and
% hence, they have an infinite number of nonzero singular values; however,
% numerically many of these are of small finite rank. For example, the
% function $f$ above is of infinite rank, but its singular values rapidly 
% decay to zero and it can be represented to machine precision by a 
% rank 5 function: 
semilogy( diag(S), '.', 'markersize', 30 )
title('The singular values of f') 

%% 
%  The singular vectors are orthogonormal with respect to the
%  inner-product. Note how we encorporate the weight function into the
%  inner-product using the "*" operator. 
rk = length(f);
weight = chebfun(@(th) sin(th),[0,pi])*ones(1,rk); 
norm( U'*(weight.*U) - eye( rk ) ) 
norm( V'*V - eye( rk ) ) 

%% Schmidt's paper from 1909
% In 1909, Schmidt wrote a wonderful paper in German entitled " ". It 
% introduced the SVD of a function, the Eckart-Young Theorem (30 years before it was 
% rediscovered by Eckart and Young []), and showed that the singular vectors 
% of a continuous function are themselves continuous. Later, work by Smithies
% paved the mathematical foundations for continuous linear algebra.  

%% 
% These early pioneers only considered the unweighted $L^2$ case of
% bivariate function on square domains. However, we translate their results 
% here (without proof) onto the sphere.   

%% 
% <latex> 
% Theorem: Let $f:S\rightarrow\mathbb{R}$ be a continuous function, then
% $f$ has a singular value decomposition in intrinsic variables, i.e.,
% 
%  $$ f(\lambda,\theta) = \sum_{k=0}^\infty \sigma_k
%  u_k(\theta)v_k(\lambda)$$ 
% 
% where each $v_k$ is a continuous function and each $u_k$ is a 
% continuous function except possibly at $\theta = 0,\pi$. 

%% 
% 
plot( U(:,3) ) 

%% The role of the SVD in Spherefun 
% Spherefun uses low rank approximations of functions to represent all
% functions and these expansions are frequently pruned and expanded as
% necessary to resolved a function to machine precision. However, 

% however, it does not compute these approximations via the SVD.
% Instead, it uses a continuous variant of Gaussian elimination that has
% been adapted to preserve symmetry on the sphere. 

%% 


%% References 
%%