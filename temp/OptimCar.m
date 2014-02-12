%% OPTIMAL PERFORMANCE OF A CAR
% Asgeir Birkisson, November 2010

%%
% (Chebfun example ode/OptimCar.m)

[d,t,N] = domain(0,2);

maxAcc = 1;

x = chebfun([0 1],d);
v = chebfun([0 1],d);
lx = chebfun([0 1],d);
lv = chebfun([0 1],d);

N.op = @(t,x,v,lx,lv) [ diff(x)-v, diff(v)-sign(1-t)*maxAcc, diff(lx),diff(lv)+lx];

N.lbc = @(x,v,lx,lv) [x,v];
N.rbc = @(x,v,lx,lv) [lx,lv];

u = N\0

plot(u),legend('Location of car','Speed','\lambda_x','\lambda_v'),
xlabel('Time'), shg