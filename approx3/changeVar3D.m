%% Change of coordinates in 3D
% Rodrigo Platte, November 2016

%%
% In this example we use mappings to compute with functions defined on 
% non-rectangular three dimensional volumes. The mapping variables must be
% defined on a rectangular domain. In other words, we will use the change
% of variables
%
% $$ x = x(u,v,w), \  y = y(u,v,w), \ z = z(u,v,w), $$
%
% where $u$, $v$, $w$ are defined as chebfun3 objects on a rectangular domain. 
%

%% Triple integrals in spherical coordinates
% Here we use spherical coordinates to compute the mass of a "ice-cream
% cone" region with variable density. The region is defined by
%

r = chebfun3(@(r,t,p) r, [0 1 0 2*pi pi/4 pi/2]);
t = chebfun3(@(r,t,p) t, [0 1 0 2*pi pi/4 pi/2]);
p = chebfun3(@(r,t,p) p, [0 1 0 2*pi pi/4 pi/2]);
x = r.*cos(t).*cos(p);
y = r.*sin(t).*cos(p);
z = r.*sin(p);

myscatter3(x,y,z), view(-53,24)

%% 
% We now define the density function and plot it

density = r.^2;
myscatter3(x,y,z,density), view(-53,24), colorbar

%%
% The mass of the solid can be found by computing the  triple integral 
% in a rectangular region. The change of variables requires us compute the
% determinant Jacobian of the transformation.

M = integral3(density.*JacDet(x,y,z));
format long
disp(M)

% BH:
M2 = integral3(density.*abs(jacobian(x,y,z)));
disp(M2)

%% 
% The exact solution is given by

disp(pi*(2-sqrt(2))/5)

%% Triple integrals in cylindrical coordinates
% In our next example we compute the center of mass of a sector of a
% cylinder with variable density.

r = chebfun3(@(r,t,z) r, [0 1 0 pi 0 1]);
t = chebfun3(@(r,t,z) t, [0 1 0 pi 0 1]);
z = chebfun3(@(r,t,z) z, [0 1 0 pi 0 1]);
x = r.*cos(t);
y = r.*sin(t); 

density = z+y;
myscatter3(x,y,z,density)
axis image, view(60,60)

%% 
% Mass:
M = integral3(density.*JacDet(x,y,z)); disp(M)

% BH:
M2 = integral3(density.*abs(jacobian(map))); disp(M2)
% or:
coord = [x; y; z];
jac = abs(jacobian(coord));
M3 = integral3(density.*jac); disp(M3)

%%
% Center of mass:
xc = integral3(x.*density.*JacDet(x,y,z))/M;
yc = integral3(y.*density.*JacDet(x,y,z))/M;
zc = integral3(z.*density.*JacDet(x,y,z))/M;
disp([xc,yc,zc])

% BH:
jac = abs(jacobian(x,y,z));
xc2 = integral3(x.*density.*jac)/M;
yc2 = integral3(y.*density.*jac)/M;
zc2 = integral3(z.*density.*jac)/M;
disp([xc2,yc2,zc2])

%% Triple integrals over the torus?
% Here is an example were we compute a triple integral over the torus

r = chebfun3(@(r,t,p) r, [0 1 0 2*pi 0 2*pi]);
t = chebfun3(@(r,t,p) t, [0 1 0 2*pi 0 2*pi]);
p = chebfun3(@(r,t,p) p, [0 1 0 2*pi 0 2*pi]);
x = (4+r.*cos(t)).*cos(p);
y = (4+r.*cos(t)).*sin(p);
z = r.*sin(t);
f = z;

%% 
myscatter3(x,y,z,f)
axis tight, axis image
view(-28,31)
disp(integral3(f.*abs(JacDet(x,y,z))))

% BH:
disp(integral3(f.*abs(jacobian(x,y,z))))



