%% Solving the diffusion equation on the surface of the unit sphere
% Alex Townsend and Grady Wright, April 2016 

%%
% (Chebfun example approx2/DiffusionEquationSphere.m)
% [Tags: #spherefun, #pdes]

%% 1. Introduction 
% Spherefun is a part of Chebfun for computing on the surface of the unit 
% sphere [1]. It presently supports about 100 commands for computing with 
% scalar- and vector-valued functions, and as has some functionality for
% solving partial differential equations.

%% 
% In this simple example we show how one can use the spherefun.Helmholtz 
% command to solve the diffusion equation on the surface of the sphere.
% We keep things simple here and just use backward Euler for the time
% integration.  Higher-order time-stepping schemes as well as more
% complicated equations, such as those involving reaction terms, should
% follow naturally from this code. 

%% 2. Heat equation on the surface of the sphere
% The heat equation defined on the surface of the sphere is given by
%
%            $$ u_t  =  \mu \nabla^2 u,  \quad     c>0, $$
%
% where $\mu$ is the coefficient of diffusion and $\nabla^2$ is the
% Laplace-Beltrami operator on the sphere. Solutions to this equation are
% unique up to an additive constant.  For an initial condition $u_0$ that
% satisfies the compatibility condition that its mean over the sphere is
% zero, this constant is typically chosen to be zero.

%% 3. Solving the diffusion equation with backward Euler
% The diffusion equation can be discretized in time using backward Euler.
% This method proceeds by discretizing the time derivative at
% $t=(k+1)\Delta t$ with the first order approximation
% $u_t\bigr|_{t=t_{k+1}} \approx(u^{k+1} - u^{k})/\Delta t$, where %t_{k} =
% k \Delta t$ and $u^{k+1}$ \& $u^{k}$ are approximations to the solution
% at time $t=t_{k+1}$ and $t_{k}$, respectively.  The Laplacian term in the
% diffusion equation is treated implicitly in time, which leads to the
% following discretization of the diffusion equation
%
%        $$  u^{k+1} - u_{k} = \mu \Delta t \nabla^2 u^{k+1},  $$ 
%
% which can be rearranged to obtain the Helmholtz equation 
% 
%       $$ \nabla^2 u^{k+1} + K^2*u^{n+1} = K^2 u^{k}, $$
% 
% where $K^2 = -\sqrt{1/(\mu \Delta t)}$. The Spherefun command
% spherefun.Helmholtz, which uses a Fourier spectral method [1], can be
% used for solving these types of equations.

%% 3. Initial condition leading to an exact solution
% As a first example, we consider an initial condition consisting of the
% following combination of spherical harmonics: 
%
% $$u_0(\lambda,\theta) = Y_{6}^0(\lambda,\theta) + \sqrt{\frac{14}{11}}
% Y_{6}^5(\lambda,\theta),$$
%
% which we can construct, and plot in Spherefun as
u0 = spherefun.sphharm(6,0) + sqrt(14/11)*spherefun.sphharm(6,5);
surf(u0), axis equal

%%
% In this case the solution of the diffusion equation at any time $t$ is
% 
% $$u(\lambda,\theta,t) = e^{-42\mu t}(Y_{6}^0(\lambda,\theta) + \sqrt{\frac{14}{11}}
% Y_{6}^5(\lambda,\theta))$$

%%
% To numerically solve the diffusion equation with this initial condition
% we use a time-step of $\Delta t = 0.01$ and integrate the solution to the
% final time $T = 1$.
dt = 0.01; 
Tfinal = 1;
t = linspace(0,Tfinal,Tfinal/dt+1); 

%%
% The code below solves the diffusion equation with $\mu = 0.02$, on a
% $32$-by-$32$ grid, specified by $n$.

mu = .02;               % Diffusion constant
n = 32;                 % Grid size
u = u0;                 % Initial condition
K = sqrt(1/dt/mu)*1i;   % Helmholtz frequency

count = 1;
for k = 2:numel(t) % each time-step
    
    % One step of backward Euler: 
    u = spherefun.Helmholtz(K^2*u, K, n, n);   % Helmholtz solve
   
    % Plot the solution at regular intervals
    if ( mod(k,25)-1 == 0 )
       subplot(2,2,count)
       plot( u ), caxis([-1 1]), title(sprintf('Time %1.2f',(k-1)*dt))
       count = count + 1;
    end    
end

%%
% Comparing the numerical solution to the exact solution at $t=1$, we find
uexact = exp(-42*mu)*u0;
norm(u-uexact)

%%
% All the error here is due to the time integrator.

%% 4. Gaussian bumps
% We conclude with a more interesting example, with an initial condition
% consisting of two Gaussian bumps and antipodal positions.
bump = @(x,y,z,x0,y0,z0,d) exp(-10*((x-x0).^2+(y-y0).^2+(z-z0).^2)); 
rng(26)
x0 = 1-2*rand(1,3); x0 = x0/norm(x0);
z0 = x0(1,3); y0 = x0(1,2); x0 = x0(1,1);
u0 = spherefun( @(x,y,z) bump(x,y,z,x0,y0,z0) + bump(x,y,z,-x0,-y0,-z0) );

%%
% To satisfy the compatibility condition of the diffusion equation, we
% adjust the initial condition so that its mean is zero.
u0 = u0 - mean2(u0);     
plot(u0), title('Initial heat profile')

%%
% We use the same time-step as the previous example, but now a value of 
% $\mu=1$ and a $100$-by-$100$ grid.
mu = 1;
n = 100;                 % Grid size
u = u0;                 % Initial condition
K = sqrt(1/dt/mu)*1i;   % Helmholtz frequency

count = 1;
for k = 2:numel(t) % each time-step
    
    % One step of backward Euler: 
    u = spherefun.Helmholtz(K^2*u, K, n, n);   % Helmholtz solve
   
    % Plot the solution at regular intervals
    if ( mod(k,25)-1 == 0 )
       subplot(2,2,count)
       contour( u ), title(sprintf('Time %1.2f',(k-1)*dt))
       count = count + 1;
    end    
end

%% 5. Poisson solver in Spherefun
% Spherefun also has a Poisson solver, which can be used to solve for 
% a steady diffusion, or heat problem; see 

help spherefun.Poisson 

%% References 
% [1] A. Townsend, H. Wilber, and G. Wright, Computing with function on
% polar and spherical geometries I. the sphere, submitted, 2016.