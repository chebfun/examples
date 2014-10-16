%% Rootrace
%
% Jared Aurentz and Anthony Austin, October 2014


%% roots of Chebyshev polynomials
%
% Here is a so called rootrace between chebfun's roots command and matlab's
% eig command.

num_degrees = 12;
degrees = 2.^(1:num_degrees)';
num_trials = flipud(degrees)/2;

% chebfun
chebfun_times = zeros(num_degrees,1);
for ii=1:num_degrees
   tic;
   for jj=1:num_trials(ii)
       %p = chebfun(randn(degrees(ii)+1,1),'coeffs');
       p = chebtech2({[], randn(degrees(ii)+1,1)});
       roots(p);
   end
   chebfun_times(ii) = toc/num_trials(ii); 
end

% matlab
matlab_times = zeros(num_degrees,1);
for ii=1:num_degrees
   oh = 0.5 * ones(degrees(ii)-1, 1);
   A = diag(oh, 1) + diag(oh, -1);
   A(end, end-1) = 1;
   e_1 = [1,zeros(1,degrees(ii)-1)]';
   tic;
   for jj=1:num_trials(ii)
       p = randn(1,degrees(ii)+1);
       eig(A-e_1*p(2:end)/p(1)/2);
   end
   matlab_times(ii) = toc/num_trials(ii); 
end

% plot times
loglog(degrees,chebfun_times,'bo',degrees,matlab_times,'r+');

% plot slopes
n = degrees/(2^6);
hold on;
loglog(degrees,n.^2,'b',degrees,n.^3,'r');
hold off;

% %% roots of Fourier polynomials
% %
% % Here is a rootrace between chebfun's roots command and matlab's
% % eig command when the functions are periodic.
% 
% % trigfun
% trigfun_times = zeros(num_degrees,1);
% for ii=1:num_degrees
%    tic;
%    for jj=1:num_trials(ii)
%        p = chebfun(randn(degrees(ii)+1,1),'coeffs','periodic');
%        roots(p);
%    end
%    trigfun_times(ii) = toc/num_trials(ii); 
% end
% 
% % matlab
% matlab_times = zeros(num_degrees,1);
% for ii=1:num_degrees
%    oh = ones(degrees(ii)-1, 1);
%    A = diag(oh, -1);
%    e_1 = [1,zeros(1,degrees(ii)-1)]';
%    tic;
%    for jj=1:num_trials(ii)
%        p = randn(1,degrees(ii)+1);
%        eig(A-e_1*p(2:end)/p(1));
%    end
%    matlab_times(ii) = toc/num_trials(ii); 
% end
% 
% % plot times
% loglog(degrees,trigfun_times,'bo',degrees,matlab_times,'r+');
% 
% % plot slopes
% n = degrees/(2^6);
% hold on;
% loglog(degrees,n.^2,'b',degrees,n.^3,'r');
% hold off;
