%% Random level hopping
% Nick Trefethen, May 2017

%%
% (Chebfun example ode-random/LevelHopping.m)

%%
% The equation $y' = -2\sin(2\pi y)$ has stable fixed points
% when $y$ is an integer.  Let us add some noise, so that we
% have
% $$ y' = -2\sin(2\pi y) + f, $$
% where $f$ is a random function.  This gives us a process
% that hops from one fixed point to another.
% We illustrate first for $t\in [0,100]$ with $\lambda = 0.4$.
rng(0), dom = [0 100]; tic
N = chebop(dom);
lambda = 0.4; f = randnfun(lambda,dom,'norm');
N.op = @(y) diff(y) + 2*sin(2*pi*y); N.lbc = 0;
LW = 'linewidth'; FS = 'fontsize';
y = N\f; plot(y,LW,0.7), grid on
xlabel('t',FS,32), ylabel('y',FS,32)

%%
% Here we cut $\lambda$ in half.
lambda = lambda/2;
f = randnfun(lambda,dom,'norm');
y = N\f; plot(y,LW,0.5), grid on
xlabel('t',FS,32), ylabel('y',FS,32)

%%
total_time_in_seconds = toc
