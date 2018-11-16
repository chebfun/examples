%% Helmholtz-Hodge decomposition of a vector field
% Nicolas Boulle and Alex Townsend, November 2018 

%% 
% (Chebfun Example ball/HelmholtzDecomposition.m)
% [Tags: #ballfun, #vector-valued, #Helmholtz-Hodge] 

%% The Helmholtz-Hodge decomposition 
% In the area of vector calculus, Helmholtz's theorem states that any
% sufficiently smooth in the ball can be expressed as a sum of a curl-free,
% a divergence-free, and an harmonic vector field [4]. In this Example, we 
% show how this decomposition is computed by Ballfun and introduce the command 
% |HelmholtzDecomposition| [3]. 

%% 
% Let $v$ be a vector field defined in the ball of radius $1$, Helmholtz's
% theorem says that we can decompose $v$ as follows: 
%
% $$ \mathbf{v} = \nabla f + \nabla \times \psi + \nabla\phi,$$
%
% where $f$ and $\phi$ are scalar-valued potential functions and $\psi$ is 
% a vector field. The first term, $\nabla f$, is a gradient field and hence
% curl-free, while the second term, \nabla \times \psi, is solenoidal. The 
% third term is an harmonic vector field (vector Laplacian of $\nabla \phi$
% is zero). From vector identities, one knows that the scalar field, $\phi$,
% is itself harmonic, i.e., $\Delta \phi = 0$.

%%
% The Helmholtz-Hodge decomposition can be made unique by imposing 
% additional constrains on $f$ and $\psi$ [4]. The standard constrains are: 
% (1) $f$ is zero on the boundary of the unit ball, (2) the normal
% component of $\psi$ on the boundary is zero, and (3) $\psi$ is 
% divergence-free. The |HelmholtzDecomposition| command in Ballfun computes 
% the decomposition under these additional constrains. 

%% 
% The Helmholtz-Hogde decomposition is an important tool in fluid dynamics as it
% is used for flow visualization (when the fluid is compressible), in CFD
% simulations (to impose that a fluid remains incompressible), and
% topological analysis. A survey of applications is available here [2]. 

%% Calculating the Helmholtz-Hodge decomposition
% To explain the procedure for computing the decomposition, we take 
% the following vector field:
v = ballfunv(@(x,y,z)cos(x.*y).*z,@(x,y,z)sin(x.*z),@(x,y,z)y.*z);
quiver( v )

%% Computing the curl-free component
% Since the divergence of a curl is zero and $\phi$ is harmonic we know that
% the divergence of $\mathbf{v}$ is the Laplacian of $f$, i.e., 
%
% $$ \nabla \cdot \mathbf{v} = \nabla \cdot \nabla f = \nabla^2 f, $$
%
% where the last equality holds because the divergence of a gradient is the
% Laplacian. Along with this, the zero Dirichlet conditions defines $f$. 

f = poisson(div(v), @(lam,th)0, 50, 50, 50);
quiver( grad( f ) ), title('Curl-free component of v')

%% 
% We confirm that this component is curl-free: 
norm( curl( grad( f ) ) ) 

%% Computing the harmonic component
% Now, one can define a vector field $v^{(1)}$ as
%
% $$v^{(1)} = v - \nabla f = \nabla \times \psi + \nabla \phi.$$ 
%
% Since we have the identity
%
% $$ \vec{n} \cdot (\nabla \times \psi)|_{\partial{B}} = 
% \psi \times \vec{n}|_{\partial B} = 0$$
%
% and want to find the harmonic function $\phi$, we solve the Laplace
% equation
%
% $$\Delta \phi = 0.$$
%
% The Neumann boundary conditions are given by
%
% $$\vec{n} \cdot \nabla \phi|_{\partial B} = \frac{\partial \phi}{\partial r}|_{\partial B}
% = \vec{n} \cdot v^{(1)}|_{\partial B}.$$
%
% Therefore, we can solve for $\psi$ in the Helmholtz-Hodge decomposition
% as follows: 
v1 = v - grad(f);
v1c = v1.comp;
Vx = coeffs3(v1c{1},50,50,50); Vy = coeffs3(v1c{2},50,50,50); Vz = coeffs3(v1c{3},50,50,50);
Vx = reshape(sum(Vx,1),50,50); Vy = reshape(sum(Vy,1),50,50); Vz = reshape(sum(Vz,1),50,50);
MsinL = trigspec.multmat(50, [0.5i;0;-0.5i] ); McosL = trigspec.multmat(50, [0.5;0;0.5] );
MsinT = trigspec.multmat(50, [0.5i;0;-0.5i] ); McosT = trigspec.multmat(50, [0.5;0;0.5] );
v1_bc = McosL*Vx*MsinT.' + MsinL*Vy*MsinT.' + Vz*McosT.';
phi = helmholtz_neumann(ballfun(@(r,lam,th)0, 'polar'), 0, v1_bc, 50, 50, 50);
quiver( grad( phi ) ), title('Harmonic component of v')

%% 
% We check the harmonicity of this component: 
norm( laplacian( grad( phi ) ) ) 

%% Computing the divergence-free component
% Let $v^{(2)}$ be the following vector field
%
% $$v^{(2)} = v^{(1)} - \nabla \phi = \nabla \times \psi.$$ 
%
% Since $v^{(2)}$ and $\psi$ are divergence-free, we can write their
% Poloidal-Toroidal decomposition [1] as
%
% $$v^{(2)} = \nabla\times\nabla\times(\mathbf{r}P_{v^{(2)}})+
% \nabla\times(\mathbf{r}T_{v^{(2)}}),$$
%
% $$\psi = \nabla\times\nabla\times(\mathbf{r}P_{\psi})+
% \nabla\times(\mathbf{r}T_{\psi}),$$
%
% where $\mathbf{r} = r\vec{r}$. Moreover, the uniqueness of the PT 
% decomposition and further vector identities lead us to the following 
% system of equations for $\psi$:
%
% $$\Delta P_{\psi} = -T_{v^{(2)}}, \quad T_{\psi} = P_{v^{(2)}},$$
%
% where $P_{\psi}$ is subjected to zero Dirichlet conditions because
% $\vec{n}\cdot\nabla\times(\mathbf{r}T_{\psi})|_{\partial B}=0$.
% Therefore, we can solve for $\psi$ in the Helmholtz-Hodge decomposition
% as follows:
v2 = v1 - grad(phi);
[Pv, Tv] = PTdecomposition(v2);
Ppsi = poisson(-Tv, @(lam,th)0, 50, 50, 50);
Tpsi = Pv;

%%
% We then recover the divergence-free vector field $\psi$ since
%
% $$\psi = \nabla\times\nabla\times(\mathbf{r}P_{\psi})+
% \nabla\times(\mathbf{r}T_{\psi}).$$
% 
psi = ballfunv.PT2ballfunv(Ppsi,Tpsi);
quiver( curl( psi ) ), title('Divergence-free component of v')

%% 
% By vector identities this component is divergence-free: 
norm( div( curl( psi ) ) ) 

%% Visualizing the decomposition
% Here is a plot of each component of the decomposition. 
subplot(1,4,1) 
quiver( v ), title('Vector field')
subplot(1,4,2)
quiver( grad(f) ), title('Curl-free')
subplot(1,4,3)
quiver( curl(psi) ), title('Divergence-free')
subplot(1,4,4)
quiver( grad(phi) ), title('Harmonic')

%% 
% As a sanity check we confirm that the decomposition has been successful:
w = grad( f ) + curl( psi ) + grad( phi );
norm( v - w ) 

%% The HelmholtzDecomposition command 
% Ballfun has a command called |HelmholtzDecomposition| that computes the 
% Helmholtz-Hodge decomposition of a vector field. Therefore, this example 
% can be replicated with the following code: 
[f, Ppsi, Tpsi, phi] = HelmholtzDecomposition( v );
psi = ballfunv.PT2ballfunv(Ppsi, Tpsi);

%% References:
%
% [1] G. Backus, Poloidal and toroidal fields in geomagnetic field modelling,
% _Reviews of Geophysics_, 24 (1986), pp. 75-109.
%
% [2] H. Bhatia, G. Norgard, V. Pascucci, and P.-T. Bremer, The 
% Helmholtz-Hodge Decomposition--A Survey, _IEEE Trans. Vis. Comput. Graphics_, 
% 19 (2013), pp. 1386-1404.
%
% [3] N. Boulle, and A. Townsend, Computing with functions in the ball, in 
% preparation.
%
% [4] Y. Tong, S. Lombeyda, A. Hirani, and M. Desbrun, Discrete Multiscale
% Vector Field Decomposition, _ACM Trans. Graphics_, 22 (2003), pp. 445-452.
