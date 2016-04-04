%% Solving the heat equation on the surface of the unit sphere
% Alex Townsend and Grady Wright, April 2016 

%% 
% Chebfun Example Spherefun folder? 

%% Introduction 
% Spherefun is a part of Chebfun for computing on the surface of the unit 
% sphere [1]. At the moment it has about 100 commands for computing with 
% scalar- and vector-valued functions. There is also some
% functionality for solving partial differential equations. 

%% 
% In this simple example we show how one can use the spherefun.Helmholtz 
% command to solve the heat equation with backward Euler.  This example 
% could easily be adapted to solve more complicated equations with high-order
% time-stepping schemes. 

%% Heat equation on the surface of the sphere
% The heat equation defined on the surface of the sphere is given by
%
%            $$ u_t  =  c*lap(u),  \quad     c>0, $$
%
% where c is the diffusion constant and lap(u) is the Laplacian of u. 
% Here, we will take the initial heat profile as a sum of Gaussian bumps:

rng(10)
init = spherefun([]); 
for bumps = 1:5
    x0 = 2*rand-1; y0 = sqrt(1-x0^2)*(2*rand-1); z0 = sqrt(1-x0^2-y0^2); 
    init = init + spherefun(@(x,y,z) exp(-30*((x-x0).^2+(y-y0).^2+(z-z0).^2)) ); 
end
init = init - mean2(init);     
plot(init), title('Initial heat profile')

%% Proceeding in time with backward Euler
% Given a time step $dt>0$, we discretize time by $0=t_0<t_1<...$ where $t_k =
% k*dt$. Let's run the heat equation for $T = 1$ seconds with $dt = 0.01$. 

dt = 0.01; 
T = 1;
t = 0:dt:T; 

%%
% For backward Euler, the laplacian term is treated implicitly and 
% $u_t$ at time $t_k$ is replaced by the first-order approximation 
% $u_t\approx(u_{k+1} - u_n)/dt$, where $u_n$ is the heat profile at time 
% $n\times dt$.  By substituting this into the heat equation we find that  
%
%        $$  u_{n+1} - u_{n} = c dt \nabla u_{n+1},  $$ 
%
% which can be rearranged to obtain the Helmholtz equation 
% 
%       $$ \nabla u_{n+1} + K^2*u_{n+1} = K^2u_{n}, $$
% 
% where K^2 = -sqrt(1/dt/c). This is the equation we solve for $u_{n+1}$ 
% at each step to proceed in time by $dt$. 

%% Heat equation solver on the surface of the unit sphere
% We can solve the heat equation in Spherefun as follows:

c = .05;           % Diffusion constant
u = init;          % Initial condition

for k = 2:numel(t) % each time-step
    
    % One step of backward Euler: 
    K = sqrt(1/dt/c)*1i;                           % Helmholtz frequency
    u = spherefun.Helmholtz(K^2*u, K, 100, 100);   % Helmholtz solve
   
    % Plot the solution at t = .25, .5, .75, and 1: 
    if ( abs( round(4*k*dt) - 4*k*dt ) < dt/100 )
       subplot(2,2, round(4*k*dt) )
       plot( u ), caxis([-1 1]), title(sprintf('Time %1.2f',k*dt))
    end
    
end

%% Poisson solver in Spherefun
% Spherefun also has a Poisson solver, see 

help spherefun.Poisson 

%% 
% This allows one to also use fully explicit or explicit-implicit
% time-stepping schemes. 

%% References 
% [1] A. Townsend, H. Wilber, and G. Wright, Computing with function on
% polar and spherical geometries I. the sphere, submitted, 2016.