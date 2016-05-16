%% Spherical harmonics
% Alex Townsend and Grady Wright, May 2016

%%
% (Chebfun example spherefun)
% [Tags: #spherefun, #eigenfunctions, #fourier]

%% 1. Introduction
% <latex>
% Spherical harmonics are the spherical analogue of trigonometric
% polynomials on the interval. The degree $\ell\geq 0$, order $m$ ($-\ell \leq m
% \leq m$) spherical harmonic is typically denoted as
% $Y_{\ell}^{m}(\lambda,\theta)$, and can be expressed in real form as [1, Sec. 14.30]:
% \begin{equation}
% Y_{\ell}^{m}(\lambda,\theta) =
% \begin{cases}
% \sqrt{2} a_{\ell}^{m} P_{\ell}^{m}(\cos\theta)\cos(m\lambda) & m > 0,\\
% a_{\ell}^{0} P_{\ell}(\cos\theta) & m=0, \\
% \sqrt{2} a_{\ell}^{|m|} P_{\ell}^{|m|}(\cos\theta)\sin(m\lambda) & m < 0,
% \end{cases}
% \end{equation}
% where $a_{\ell}^{k}$, $0\leq k \leq \ell$, is a normalization factor
% and $P_{\ell}^{k}$, $0\leq k \leq \ell$, is the degree $\ell$ and
% order $k$ associated Legendre functions [1, Ch. 14].
% Here, we have used the following spherical coordinate parameterization 
% for a point on the unit sphere $(x,y,z)$:
% % $$ x = \cos\lambda\sin\theta,\; y = \sin\lambda\sin\theta,\; z =
% \cos\theta,$$
% where $-\pi \leq lambda \leq \pi$ and $0 \leq \theta \leq \pi$.
% </latex>

%%
% <latex>
% Spherical harmonics can be derived by solving the eigenvalue problem
% for the surface Laplace (Laplace-Beltrami) operator on the sphere; for 
% an alternative derivation see [2, Ch.\ 2].
% This operator can be expressed in the spherical coordinates defined above
% as
% $$ \Delta = \frac{\partial^2}{\partial \theta^2}
% - \frac{\cos\theta}{\sin\theta}\frac{\partial}{\partial\theta} +
% \frac{1}{\sin^2\theta}\frac{\partial^2}{\partial\lambda^2}. $$
% Any degree $\ell\geq 0$ spherical harmonic has the property that 
% $\Delta Y_{\ell}^{m} = -\ell(\ell + 1) Y_{\ell}^{m}$, 
% $-m\leq \ell \leq m$.  So, there are 2\ell+1 eigenfunctions of $\Delta$ 
% with eigenvalue $-\ell(\ell+1)$. Spherical harmonics can also be
% expressed in Cartesian form as polynomials of $x$, $y$, and $z$ [2, Ch.
% 2]. When viewed in this way, one finds that these polynomials all 
% satisfy Laplace's equation in $\mathbb{R}^3$, i.e., they are _harmonic.
% This is where the name spherical harmonics originates and was first used
% by Thompson (Lord Kelvin) and Tait [3, Appendix B].
% </latex>

%%
% <latex>
% By choosing the normalization factors $a_{\ell}^{k}$ in (1) as 
% $$ a_{\ell}^{k} = \sqrt{\frac{(2\ell+1)(\ell-k)!}{4\pi(\ell+m)!}}, $$
% the set of spherical harmonics $\{Y_{\ell}^{m}\}$, $\ell=0,1,\ldots$, 
% $m = -\ell,\ldots,\ell$, is orthonormal, i.e., 
% $$ \int_{-\pi}^{\pi} \int_0^{\pi}
% Y_{ell}^{m}(\lambda,\theta)Y_{ell'}^{m'}(\lambda,\theta)
% \sin\theta\,d\theta d\lambda .$$
% Furthermore, it can be shown that they form a complete orthonormal basis
% for the set of all $L^{2}$ integrable functions on the sphere, 
% $L^{2}(\mathbb{S})$, [2, Sec. 2.8]. For any $f\inL^{2}(\mathbb{S})$, we
% have 
% \begin{equation}
% f(\lambda,\theta) = \sum_{\ell=0}^{\infty\sum_{m=-\ell}^\ell c_{\ell,m} Y_{\ell}^{m}(\lambda,\theta),
% \label{eq:sphHarmEx}
% \end{equation}
% where 
% \begin{equation}
% c_{\ell,m} = \int_{-\pi}^{\pi} \int_0^{\pi}
% f(\lambda,\theta)Y_{ell}^{m}(\lambda,\theta)\,d\theta d\lambda,
% \end{equation}
% and equality in (\ref{eq:sphHarmEx}) is understood in the mean-square 
% sense. Truncating the outer sum of (\ref{eq:sphHarmEx}) to $N$ gives 
% the best degree in $N$ approximation of $f$ in the $L^2$ norm on the 
% sphere amongst all harmonic polynomials in $\mathbb{R}^3$ of degree 
% $N$ restricted to the sphere [2, Ch. 4].  A spherical harmonic expansion
% gives essentially uniform resolution of a function over the sphere.
% </latex>

%% 2. Spherical harmonics in Spherefun
% While spherical harmonics have many properties that make them 
% mathematically appealing to represent functions on the sphere, the
% technology behind Spherefun does not rely on them.  Instead, it combines
% the double Fourier Sphere sphere method with a low rank 
% approximation method (based on structure preserving Gaussian elimination) 
% for approximating functions on the sphere to approximately machine 
% precision [4]. This latter method allows for fast, highly adaptable
% discretizations based on the FFT, with no setup cost.  While fast
% spherical harmonic transforms are available, the precomputation cost for
% these algorithms is too high to use to meet the requirements of 
% Spherefun.  

%%
% Nevertheless, given the dominance of spherical harmonics in many
% applications, Spherefun provides many commands that can be used for
% computing with spherical harmonics.  First amongst these is the |sphharm|
% function.  For example, $Y_{17}^{13}$ can be constructed and plotted in
% Spherefun using
Y17 = spherefun.sphharm(17,13);
plot(Y17)

%%
% We can verify that this function is an eigenfunction of the surface 
% Laplacian with eigenvalue $-\ell(\ell+1) = 17(18)$:
norm(laplacian(Y17)-(-17*18)*Y17)

%%
% We can also verify the orthogonality of spherical harmonics on the sphere
Y13 = spherefun.sphharm(13,7);
sum2(Y13.*Y17)
sum2(Y13.*Y13)
sum2(Y17.*Y17)

%% 3. Computing spherical harmonic coefficients
% While we have implemented a fast spherical harmonic transform in
% Spherefun, one can use the |sum2| command to compute the spherical
% harmonic coefficients. Classically, for a band-limited function $f$ with 
% a bandwidth of $N$ this costs a totally $O(N^4)$ operations; however, 
% in Spherefun if $f$ is of low rank then this is reduced to $O(N^3)$ 
% operations. Here, we take a Gaussian bump on the sphere and then complete 
% some of its low spherical harmonic modes. 

% Example: Gaussian
rng(10)
x0 = 2*rand-1; y0 = sqrt(1-x0^2)*(2*rand-1); z0 = sqrt(1-x0^2-y0^2);
sig = 0.05;
f = spherefun(@(x,y,z) exp(-((x-x0).^2+(y-y0).^2+(z-z0).^2))/sig^2 )
plot(f), title('A Gaussian bump on the sphere',FS,16)

%%
% Spherical harmonic coefficients are computed as
k = 1;
for ell = 0:10
    for m = -ell:ell
        Y = spherefun.sphharm(ell,m);
        coeffs(k) = sum2(f.*Y);
        b(k) = Y(x0,y0,z0);
        k = k + 1;
    end
end
semilogy(coeffs./b,'x-')
title('The spherical harmonic coefficients of a Gaussian bump',FS,16)
set(gca,FS,16)

%% 
% Why do we get the strange step effect??   

%% 
% Exact values seems hard... 

%% 4. Platonic solids
% Certain low order spherical harmonics can be combined so that they have
% the same rotational symmetries of certain platonic solids.  These
% combinations of spherical harmonics play a key role in the linear
% stability analysis of partial differential equations in spherical
% geometries; see, for example, the work of Busse on convection in
% spherical shells [5].  The combination with tetrahedral symmetry is 
% given by
Y = spherefun.sphharm(3,2);
plot(Y), colormap(jet(2))

%%
% The combination with octahedral symmetry is given by
Y = spherefun.sphharm(4,0) + sqrt(5/7)*spherefun.sphharm(4,4);
plot(Y)

%%
% And the combination with icosahedral symmetry is given by
Y = spherefun.sphharm(6,0) + sqrt(14/11)*spherefun.sphharm(6,5);
plot(Y)

%% References
% [1] F. W. Olver, D. W. Lozier, R. F. Boisver, and C. W. Clark, _NIST
% Handbook of Mathematical Functions_, Cambridge University Press, 2010.
%%
% [2] K. Atkinson and W. Han, _Spherical Harmonics and Approximations on the Unit Sphere: An
% Introduction_, Lecture Notes in Mathematics, Springer, 2012.
%%
% [3] W. Thomson and P. G. Tait. Treatise on Natural Philosophy. 2nd ed.
% Vol. 1, Cambridge: At the University press, 1888.
%%
% [4] A. Townsend, H. Wilber, and G. B. Wright, Computing with function in
% polar and spherical geometries I. The sphere, to appear in 
% _SIAM J. Sci. Comp._, 2016 
%%
% [5] F. H. Busse, Patterns of convection in spherical shells. _J. Fluid
% Mech._, 72, 67-85, 1975.