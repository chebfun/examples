%% Odd and even best approximations
% Mohsin Javed and Nick Trefethen, March 2015

%%
% (Chebfun example approx/OddEven.m)
% [Tags: #REMEZ, #bestapproximation, #oddeven]

%% 
% To find the best (minimax) approximation of a function
% $f$, can you find best approximations to the even part and
% the odd part, and add them together?

%%
% Such additivity would certainly apply for a linear approximation
% process such as interpolation in Chebyshev points.
% Since best approximation is nonlinear, however, one would
% expect that the additivity would fail. This is indeed the case,
% as we can easily show with an example.
% Here is a Gaussian defined on $[-1,1]$ and its
% best approximant of degree $0$:
f = chebfun(@(x) exp(-150*(x-.5).^2));
[p,err] = remez(f,0);
LW = 'linewidth'; lw = 1.6; FS = 'fontsize'; fs = 10;
CO = 'color'; green = [0 .7 0]; ax = [-1 1 -1.2 1.2];
plot(f,'b',p,'r',LW,lw), axis(ax)
grid on, title('f and its best approximation',FS,fs)

%%
% The error looks like this:
plot(p-f,CO,green,LW,lw), grid on
title(['error = ' num2str(err)],FS,fs), axis(ax)

%%
% Here is the even part of $f$ and its best approximant:
x = chebfun('x');
feven = (f(x) + f(-x))/2; [peven,erreven] = remez(feven,0);
subplot(2,2,1), plot(feven,'b',peven,'r',LW,lw), axis(ax)
grid on, title('Approximation of the even part',FS,fs)
subplot(2,2,2), plot(peven-feven,CO,green,LW,lw), grid on
title(['error = ' num2str(erreven)],FS,fs), axis(ax)

%%
% Here is the odd part of $f$ and its best approximant,
% namely the zero function.
fodd = (f(x) - f(-x))/2; [podd,errodd] = remez(fodd,0);
subplot(2,2,1), plot(fodd,'b',podd,'r',LW,lw), axis(ax)
grid on, title('Approximation of the odd part',FS,fs)
subplot(2,2,2), plot(podd-fodd,CO,green,LW,lw)
grid on, title(['error = ' num2str(errodd)],FS,fs), axis(ax)

%%
% Now, if we add up the even approximation and the odd
% approximation, how does the combination do?
% We see that the error is greater than before.
psum = peven + podd;
errsum = norm(f-psum,inf);
subplot(2,2,1), plot(f,'b',psum,'r',LW,lw), axis(ax)
grid on, title('combined',FS,fs)
subplot(2,2,2), plot(psum-f,CO,green,LW,lw), grid on
title(['error = ' num2str(errsum)],FS,fs), axis(ax)

%%
% Here is a second example, but with approximations of degree
% 1.  To ensure there are enough oscillation points to make
% the best approximations elegant, we upgrade
% our camel from dromedary to bactrian:
f = chebfun(@(x) exp(-300*(x-.25).^2) + exp(-300*(x-.75).^2));
[p,err] = remez(f,1);
subplot(2,2,1), plot(f,'b',p,'r',LW,lw)
grid on, title('f and its best approximation',FS,fs), axis(ax)
subplot(2,2,2), plot(p-f,CO,green,LW,lw), grid on
title(['error = ' num2str(err)],FS,fs), axis(ax)

%%
% The even part and its best approximant now look like this,
feven = (f(x) + f(-x))/2; [peven,erreven] = remez(feven,1);
subplot(2,2,1), plot(feven,'b',peven,'r',LW,lw)
grid on, title('Approximation of the even part',FS,fs), axis(ax)
subplot(2,2,2), plot(peven-feven,CO,green,LW,lw), grid on
title(['error = ' num2str(erreven)],FS,fs), axis(ax)

%%
% and the odd part and its best approximation look like this,
fodd = (f(x) - f(-x))/2; [podd,errodd] = remez(fodd,1);
subplot(2,2,1), plot(fodd,'b',podd,'r',LW,lw)
grid on, title('Approximation of the odd part',FS,fs), axis(ax)
subplot(2,2,2), plot(podd-fodd,CO,green,LW,lw), grid on
title(['error = ' num2str(errodd)],FS,fs), axis(ax)

%%
% Again, the sum of the two is not as good an approximation.
psum = peven + podd; errsum = norm(f-psum,inf);
subplot(2,2,1), plot(f,'b',psum,'r',LW,lw)
grid on, title('combined',FS,fs), axis(ax)
subplot(2,2,2), plot(psum-f,CO,green,LW,lw), grid on
title(['error = ' num2str(errsum)],FS,fs), axis(ax)

%%
% The reader may enjoy exploring other functions $f$ with
% the code given here.
% If $f$ is neither even nor odd, you will probably find as we have
% that some optimality is lost when the even and odd parts are separated.
