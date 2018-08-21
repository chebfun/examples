%% Helmholtz-Hodge decomposition of a vector field
% Nicolas Boull» and Alex Townsend, August 2018 

%% 
% (Chebfun Example ball/HelmholtzDecomposition.m)
% [Tags: #ballfun, #vector-valued, #Helmholtz-Hodge] 

%% 1. The Helmholtz-Hodge decomposition 
% According to the Helmholtz-Hodge theorem [1], any smooth vector field in 
% the ball can be uniquely decomposed into a sum of a divergence-free component,
% a curl-free component and a harmonic component.
% In this Example, we show how this decomposition can be computed by
% Ballfun and introduce the command |HelmholtzDecomposition| [2]. 

%% 
% Let $v$ be any vector field in the ball, we can write it as the sum
% $$ \mathbf{v} = \nabla f + \nabla \times \psi + \nabla \phi, $$
% where $f$ and $\phi$ are scalar-valued potential functions and $\psi$ is 
% a vector field. Moreover, $f|_{\partial B} = 0$,
% $\psi\times\vec{n}|_{\partial{B}} = 0$, $\psi$ is divergence-free and
% $\phi$ is harmonic is $\Delta \phi = 0$. Here $\nabla$ is the gradient 
% on the ball, $\nabla \times \psi$ stands for the divergence of $\psi$ and
% $\vec{n}$ is the unit normal vector.

%% 
% To illustrate the decomposition, we take the following 
% vector field:
S = [40,40,40];
vx = ballfun(@(x,y,z)cos(x.*y).*z,'cart',S);
vy = ballfun(@(x,y,z)sin(x.*z),'cart',S);
vz = ballfun(@(x,y,z)y.*z,'cart',S);
v = ballfunv(vx,vy,vz);
quiver( v ), view([-36 8])

%% 2. Computing the curl-free component
% Since the divergence of a curl is zero and $\phi$ is harmonic we know that 
% $$ \nabla \cdot \mathbf{c} = \nabla \cdot \nabla \psi = \nabla^2 \phi, $$
% where the last equality holds because the divergence of a gradient is the
% Laplacian. Therefore, we impose homogeneous boundary conditions on $f$ and
% solve the equation as follows: 
f = helmholtz(div(v),0,zeros(40,40));
quiver( grad( f ) ), hold on,
title('Curl-free component of v')
view([-36 8]), hold off

%% 
% We confirm that this component is curl-free: 
norm( curl( grad( f ) ) ) 

%% 3. Computing the harmonic component
% We define a vector field $v^{(1)}$ as
% $$v^{(1)} = v - \nabla f = \nabla \times \psi + \nabla \phi.$$ 
% Since we have the identity
% $$ \vec{n} \cdot (\nabla \times \psi)|_{\partial{B}} = 
% \psi \times \vec{n}|_{\partial B} = 0$$
% and want to find the harmonic component $\phi$, we solve the Laplace
% equation
% $$\Delta \phi = 0.$$
% The Neumann boundary condition is given by
% $$\vec{n} \cdot \nabla \phi|_{\partial B} = \frac{\partial \phi}{\partial r}|_{\partial B}
% = \vec{n} \cdot v^{(1)}|_{\partial B}.$$
% Therefore, we can solve for $\psi$ in the Helmholtz-Hodge decomposition
% as follows: 
v_1 = v - grad(f);
v_1_Boundary = ComputeNormalBoundary(v_1);
phi = helmholtz_neumann(ballfun(zeros(40,40,40)),0,v_1_Boundary);
quiver( grad( phi ) ), hold on,
title('harmonic component of v')
view([-36 8]), hold off

%% 
% We check the harmonicity of this component: 
norm( laplace( grad( phi ) ) ) 

%% 3. Computing the divergence-free component
% Let $v^{(2)}$ be the following vector field
% $$v^{(2)} = v^{(1)} - \nabla \phi = \nabla \times \psi.$$ 
% Since $v^{(2)}$ and $\psi$ are divergence-free, we can write their
% Poloidal-Toroidal decomposition 
% $$v^{(2)} = \nabla\times\nabla\times(\mathbb{r}P_{v^{(2)}})+
% \nabla\times(\mathbb{r}T_{v^{(2)}}),$$
% $$\psi = \nabla\times\nabla\times(\mathbb{r}P_{\psi})+
% \nabla\times(\mathbb{r}T_{\psi}).$$
% Moreover, the uniqueness of the PT decomposition and the vector
% identities lead us to the following system of equations for $\psi$:
% $$\Delta P_{\psi} = -T_{v^{(2)}}, \quad T_{\psi} = P_{v^{(2)}},$$
% where the first equation is subjected to the Dirichlet boundary condition
% $P_{\psi}|_{\partial B} = 0$ as
% $\vec{n}\cdot\nabla\times(\mathbb{r}T_{\psi})|_{\partial B}=0.$
% Therefore, we can solve for $\psi$ in the Helmholtz-Hodge decomposition
% as follows:
v_2 = v_1 - grad(phi);
[Pv,Tv] = PTdecomposition(v_2);
Ppsi = helmholtz(-Tv,0,zeros(40,40));
Tpsi = Pv;
% We then recover the divergence-free vector field $\psi$ since
% $$\psi = \nabla\times\nabla\times(\mathbb{r}P_{\psi})+
% \nabla\times(\mathbb{r}T_{\psi}).$$
psi = ballfunv.PT2ballfunv(Ppsi,Tpsi);
quiver( curl( psi ) ), hold on,
title('divergence-free component of v')
view([-36 8]), hold off

%% 
% By vector identities this component is divergence-free: 
norm( div( curl( psi ) ) ) 

%% 4. Plotting the decomposition
% Here is a plot of the decomposition. 
subplot(1,4,1) 
quiver( v ), title('vector field'), view([-36 8])
subplot(1,4,2)
quiver( grad(f) ), title('Curl-free'), view([-36 8])
subplot(1,4,3)
quiver( curl(psi) ), title('Divergence-free'), view([-36 8])
subplot(1,4,4)
quiver( grad(phi) ), title('Harmonic'), view([-36 8])

%% 
% As a sanity check we confirm that the decomposition has been successful:
w = grad( f ) + curl( psi ) + grad(phi);
norm( v - w ) 

%% 5. The HelmholtzDecomposition command 
% Ballfun has a command called |HelmholtzDecomposition| that
% computes the Helmholtz-Hodge decomposition of a vector field.
% Therefore, this example can be replicated with the following code: 
[f,Ppsi,Tpsi,phi] = HelmholtzDecomposition(v);
psi = ballfunv.PT2ballfunv(Ppsi,Tpsi);

%% 6. References
%
% [1] Y. Tong, S. Lombeyda, A. Hirani, and M. Desbrun, Discrete Multiscale
% Vector Field Decomposition, _ACM Trans. Graphics_, 22 (2003), pp. 445-452.
% 
% [2] N. Boull», J. Slomka, and A. Townsend, Optimal complexity Navier-Stokes 
% simulations in the ball, in preparation.
