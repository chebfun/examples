%% Extrema of a complicated function
% Nick Trefethen, September 2010

%%
% (Chebfun example opt/ExtremeExtrema.m)
% [Tags: #optimization, #localextrema, #MIN, #MAX]

%%
% Here is the function $\cos(x)\sin(e^x)$ on the interval $[0,6]$:
tic, x = chebfun('x',[0 6]);
f = cos(x).*sin(exp(x));
LW = 'linewidth';
length(f)
plot(f,LW,1,'color',[0 .7 0])
FS = 'fontsize'; MS = 'markersize';
title('A complicated function',FS,12)

%%
% Here's its absolute value:
g = abs(f);
ax = [0 6 0 1];
plot(g,'m',LW,1), axis(ax)
title('Absolute value',FS,12)

%%
% Here's the minimum of that function and $x/8$:
h = min(g,x/8);
plot(h,LW,1), axis(ax)
title('Minimum with x/8',FS,12)

%%
% We can find the maximum over the interval $[0,5]$ like this:
[maxval,maxpos] = max(h{0,5})
hold on, plot(maxpos,maxval,'.r',MS,30)
title('Global maximum',FS,12)

%%
% Let's add all the local maxima to the plot:
[val,pos] = max(h,'local');
plot(pos,val,'.k',MS,12)
plot(maxpos,maxval,'.r',MS,30)
title('Local maxima',FS,12)

%%
% These computations showcase the fact that Chebfun optimization is global --
% whether in the sense of finding a global extremum, or in the sense of
% globally finding all the local extrema.

%%
% They also showcase the treatment of discontinuities. To find extrema,
% Chebfun examines zeros of the derivative. In some cases those are points
% where the derivative is continuous and passes through zero. In others,
% like black dot near $x=1$
% in the plot above, they are points where the
% derivative jumps from positive to negative or vice versa.

%%
% Here is the time for this whole sequence of computations:
Total_time = toc
