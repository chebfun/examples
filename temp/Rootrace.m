%% Rootrace
%
% Jared Aurentz, 16 October 2014

%% Introduction
% The core of Chebfun rootfinding is recursive subdivision. This technique
% reduces the complexity of Matlab's |roots| (|eig|) command from $O(n^3)$ 
% to $O(n^2)$. An important question is, _How big are the constants?_ To 
% answer this question the roots of a sequence of polynomials with random 
% coefficents and increasing are computed using both Chebfun and Matlab.
% The average computation times are then compared graphically.

%% Roots of Chebyshev polynomials
%
% In the first set of examples roots of polynomials expressed in a 
% Chebyshev basis with normally distributed real coefficients are computed
% using two different methods. The first method constructs a |chebfun|
% object and calls the overloaded Chebfun |roots|. The second method
% constructs the corresponding colleague matrix and calls Matlab's |eig|
% command.

clear all
close all

num_degrees = 12;
degrees = 2.^(1:num_degrees)';
num_trials = flipud(degrees)/2;

chebfun_times = zeros(num_degrees,1);
for ii=1:num_degrees
   for jj=1:num_trials(ii)
       p = chebfun(randn(degrees(ii)+1,1),'coeffs');
       tic;
       roots(p);
       chebfun_times(ii) = chebfun_times(ii) + toc;
   end
   chebfun_times(ii) = chebfun_times(ii)/num_trials(ii); 
end

matlab_times = zeros(num_degrees,1);
for ii=1:num_degrees
   oh = 0.5 * ones(degrees(ii)-1, 1);
   A = diag(oh, 1) + diag(oh, -1);
   A(end, end-1) = 1;
   e_1 = [1,zeros(1,degrees(ii)-1)]';
   for jj=1:num_trials(ii)
       p = randn(1,degrees(ii)+1);
       tic;
       eig(A-e_1*p(2:end)/p(1)/2);
       matlab_times(ii) = matlab_times(ii) + toc;
   end
   matlab_times(ii) = matlab_times(ii)/num_trials(ii); 
end

MS = 'markersize'; FS = 'fontsize'; LW = 'linewidth';
loglog(degrees,chebfun_times,'bo',degrees,matlab_times,'r+',MS,10,LW,2);
title('Chebfun v. Matlab',FS,16);
xlabel('degree',FS,16);
ylabel('average time (sec)',FS,16);
legend('Chebfun','Matlab','location','NorthWest');

%%
% The speed up at the highest degree is pretty good.
degree = degrees(end)
speedup = matlab_times(end)/chebfun_times(end)

%% 
% For the higher degrees the asymptotic behavior becomes apparent. For the
% lower degrees the overhead of calling Chebfun's |roots| command dwarfs the
% cost of the actual root finding algorithm. We can strip off some of this
% overhead by using the lower level class |chebtech2|.

chebtech_times = zeros(num_degrees,1);
for ii=1:num_degrees
   for jj=1:num_trials(ii)
       p = chebtech2({[], randn(degrees(ii)+1,1)});
       tic;
       roots(p);
       chebtech_times(ii) = chebtech_times(ii) + toc;
   end
   chebtech_times(ii) = chebtech_times(ii)/num_trials(ii); 
end

loglog(degrees,chebfun_times,'bo',degrees,chebtech_times,'kx',degrees,matlab_times,'r+',MS,10,LW,2);
title('Chebfun v. Chebtech v. Matlab',FS,16);
xlabel('degree',FS,16);
ylabel('average time (sec)',FS,16);
legend('Chebfun','Chebtech','Matlab','location','NorthWest');

%%
% If you look closely you can see a small jump in the computation times for
% |chebtech2| at degree 50. This is where the recursive subdivision begins to
% take effect.

%% Roots of Fourier polynomials
%
% Here the same experiment is performed for trigonometric polynomials.

trigfun_times = zeros(num_degrees,1);
for ii=1:num_degrees
   for jj=1:num_trials(ii)
       p = chebfun(randn(degrees(ii)+1,1),'coeffs','periodic');
       tic;
       roots(p);
       trigfun_times(ii) = trigfun_times(ii) + toc;
   end
   trigfun_times(ii) = trigfun_times(ii)/num_trials(ii); 
end

trigtech_times = zeros(num_degrees,1);
for ii=1:num_degrees
   for jj=1:num_trials(ii)
       p = trigtech({[], randn(degrees(ii)+1,1)});
       tic;
       roots(p);
       trigtech_times(ii) = trigtech_times(ii) + toc;
   end
   trigtech_times(ii) = trigtech_times(ii)/num_trials(ii); 
end

matlab_times = zeros(num_degrees,1);
for ii=1:num_degrees
   for jj=1:num_trials(ii)
       p = randn(1,degrees(ii)+1);
       tic;
       roots(p);
       matlab_times(ii) = matlab_times(ii) + toc;
   end
   matlab_times(ii) = matlab_times(ii)/num_trials(ii); 
end

loglog(degrees,trigfun_times,'bo',degrees,trigtech_times,'kx',degrees,matlab_times,'r+',MS,10,LW,2);
title('Trigfun v. Trigtech v. Matlab',FS,16);
xlabel('degree',FS,16);
ylabel('average time (sec)',FS,16);
legend('Trigfun','Trigtech','Matlab','location','NorthWest');

%%
% The speed up at the highest degree for |trigfun|.
degree = degrees(end)
speedup = matlab_times(end)/trigfun_times(end)
