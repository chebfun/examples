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
% \label{eq:sphCoeffs}
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
% spherical harmonic transforms are available [5], the precomputation cost
% for these algorithms is too high to meet the requirements of Spherefun.

%%
% Nevertheless, given the dominance of spherical harmonics in many
% applications, Spherefun provides many commands that can be used for
% computing with spherical harmonics.  First amongst these is the |sphharm|
% function, which is used to construct a spherical harmonic of a given
% degree and order. For example, $Y_{17}^{13}$ can be constructed and 
% plotted as follows:
Y17 = spherefun.sphharm(17,13);
plot(Y17)

%%
% We can verify that this function is an eigenfunction of the surface 
% Laplacian with eigenvalue $-\ell(\ell+1) = 17(18)$:
norm(laplacian(Y17)-(-17*18)*Y17)

%%
% We can also verify the orthogonality of spherical harmonics on the sphere
% using |sum2|, which computes the integral of a function over the sphere: 
Y13 = spherefun.sphharm(13,7);
sum2(Y13.*Y17)
sum2(Y13.*Y13)
sum2(Y17.*Y17)

% %% 
% % The spherical harmonics of degree $\ell \leq 2$ 
% N = 2;
% for ell = 0:N
%     for m = -ell:ell
%         Y = spherefun.sphharm(ell,m);
%         subplot(N+1,2*N+1,ell*(2*N+1)+(N+1)+m), plot(Y), axis off
%     end
% end

%% 
% Here is a plot of the real spherical harmonics $Y_{\ell}^{m}$, 
% $\ell=0,\ldots,4$ and $0\leq m \leq \ell$, with black contour lines
% marking the lines where they are zero
N = 4;
for ell = 0:N
    for m = 0:ell
        Y = spherefun.sphharm(ell,m);
        subplot(N+1,N+1,ell*(N+1)+m+1), plot(Y), hold on 
        contour(Y,[0 0],'k-'), axis off, hold off
    end
end

%%
% The negative order real spherical harmonics are similar to the positve 
% order ones and only differ by a rotation about the north and 
% south poles.

%% 3. Computing spherical harmonic coefficients
% The computational cost of computing all the spherical harmonic
% coefficients up to degree $N$ of a function directly using an
% approximation of \ref{eq:sphCoeffs} scales like $O(N^4)$.  If $f$ is low
% rank, then the coefficients can be obtained in $O(N^3)$ operations using
% a fast multiplication algorithm and the |sum2| command in Spherefun. As
% noted above, fast $O(N^2\log N)$ algorithms are available for this task
% [5], but these are not yet available in Spherefun.

%%
% As an example, consider the restriction of a Gaussian function, 
% $$ f(x,y,z) = \exp\left(-\frac{(x-x_0)^2 + (y-y_0)^2 + (z-z_0)^2}{\sigma^2}\right),$$
% to the surface of the sphere.  Here we center the Gaussian at a 
% random point $(x_0,y_0,z_0)$ on the sphere:
rng(10)
x0 = 2*rand-1; y0 = sqrt(1-x0^2)*(2*rand-1); z0 = sqrt(1-x0^2-y0^2);
sig = 0.5;
f = spherefun(@(x,y,z) exp(-((x-x0).^2+(y-y0).^2+(z-z0).^2)/sig^2) )
plot(f), title('A Gaussian bump on the sphere')

%%
% The numerical rank of this function is only 14, so we consider it a low
% rank function. The spherical harmonic coefficients of $f$ up to degree 12
% can be computed as follows:
N = 12; k = 1;
coeffs = zeros((N+1)^2,3); shifts = coeffs(:,1);
for ell = 0:N
    for m = -ell:ell
        Y = spherefun.sphharm(ell,m);
        coeffs(k,1) = sum2(f.*Y);
        coeffs(k,2:3) = [ell m];
        shifts(k) = Y(x0,y0,z0);
        k = k + 1;
    end
end
%%
% Since the Gaussian is analytic, the spherical harmonic coefficients decay
% exponentially fast with increasing $\ell$ [2], as the following figure
% illustrates:
stem3(coeffs(:,2),coeffs(:,3),abs(coeffs(:,1)),'filled'), ylim([-N N])
set(gca,'ZScale','log'), set(gca,'Xdir','reverse'), view([-13 18])
xlabel('$\ell$','Interpreter','Latex'), ylabel('m'), zlabel('|coeffs|')

%%
% Here is the least squares projection of degree 6 of the Gaussian 
% function given above:
fproj = spherefun([]);
k = 1;
for ell = 0:6
    for m = -ell:ell
        fproj = fproj + coeffs(k,1)*spherefun.sphharm(ell,m);
        k = k + 1;
    end
end
plot(fproj), title('Degree 6 spherical harmonic projection of the Gaussian');

%%
% The error in this approximation looks as follows
plot(f-fproj), colorbar, title('Error in the spherical harmonic projection')
norm(f-fproj)

%% Postive definite zonal kernels
% Normalize spherical harmonic coefficients of the Gaussian
ncoeffs = coeffs(:,1)./shifts;
stem3(coeffs(:,2),coeffs(:,3),abs(ncoeffs),'filled')
set(gca,'ZScale','log'), set(gca,'Xdir','reverse'), view([-13 18])
xlabel('$\ell$','Interpreter','Latex'), ylabel('m'), zlabel('|coeffs|')

%% 
% Why do we get the strange step effect??  Funk-Hecke formula for zonal
% kernel

%% 
% The exact values at each step are [6]
cexact = zeros((N+1)^2,1);
k = 1;
for ell = 0:N
    for m = -ell:ell
        cexact(k) = 2*sqrt(pi)^3*sig*exp(-2/sig^2).*besseli(ell+1/2,2/sig^2);
        k = k + 1;
    end
 end

%%
% Check
norm(ncoeffs(:,1)-cexact)

%% 4. Platonic solids
% Certain low order spherical harmonics can be combined so that they have
% the same rotational symmetries of certain platonic solids.  These
% combinations of spherical harmonics play a key role in the linear
% stability analysis of partial differential equations in spherical
% geometries; see, for example, the work of Busse on convection in
% spherical shells [7].  The combination with tetrahedral symmetry is 
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
% [5] M. Tygert, Fast algorithms for spherical harmonic expansions, III, _
% J. Comput. Phys._, 229 , 6181?6192, 2010.
%%
% [6] Hubbert, S. and Baxter, B., _Radial basis functions for the sphere_, 
% Progress in Multivariate Approximation, Volume 137 of the International
% Series of Numerical Mathematics, Birkhauser, 33-47, 2001.
%%
% [7] F. H. Busse, Patterns of convection in spherical shells. _J. Fluid
% Mech._, 72, 67-85, 1975.