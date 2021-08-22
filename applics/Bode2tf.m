%% The AAA algorithm for system identification
% Stefano Costa, August 2021

%%
% (Chebfun example applics/Bode2tf.m)
% [Tags: #AAA]

%%
% The AAA algorithm provides a natural way to identify LTI (linear time-invariant)
% system parameters such as
% poles, zeros and DC gain from Bode plots. For example, consider the 4th order system
% $$ G(z) = 2\frac{(1+105z)(1+1.4/0.05z+1/0.05^2 z^2)}{(1+100z)(1+1.4/0.04z+1/0.04^2 z^2)(1+z)} $$
Nc = 2*conv([105 1],[1/0.05^2 1.4/0.05 1]);
Dc = conv(conv([100 1],[1/0.04^2 1.4/0.04 1]),[1 1]);
N = @(z) Nc*z.^[3:-1:0]';
D = @(z) Dc*z.^[4:-1:0]';
G = @(z) N(z(:))./D(z(:));

%%
% for which we have the following:
pol = roots(Dc), zer = roots(Nc), DCgain = abs(G(0))

%%
% Let's sample some values in the frequency range $10^{-4}\leq \omega \leq 10^2$, and
% draw the Bode plots of magnitude and phase:
w = logspace(-4,2,3000);
mag = abs(G(i*w)); ph = -angle(G(i*w));
LW = 'linewidth'; LO = 'location'; SW = 'southwest';
subplot(211), semilogx(w,20*log10(mag),'b-',LW,1.5), grid on, title('Magnitude (dB)')
subplot(212), semilogx(w,ph*180/pi,'b-',LW,1.5), grid on, title('Phase (degrees)')

%%
% Given magnitude and phase, an approximation $H(s)$ for $G(s)$ is readily obtained by
% AAA approximation of the complex signal. Samples are mirrored in order to enforce
% symmetry.
wA = [-fliplr(w) w]; magA = [fliplr(mag) mag]; phA = [-fliplr(ph) ph];
GA = magA.*exp(i*phA);
[H,polA,resA,zerA] = aaa(GA,i*wA); polA, zerA, DCgainA = abs(H(0))
subplot(211), hold on, semilogx(w,20*log10(H(i*w)),'k--',LW,1.5)
subplot(212), hold on, semilogx(w,angle(H(i*w))*180/pi,'k--',LW,1.5)

%%
% Also, $H(s)$ features negligible errors in initial data:
err_mag = norm(mag-abs(H(i*w)),inf)
err_ph = norm(ph-angle(H(i*w)),inf)

%% 
% The following means of recomputing poles and zeros will play a key role in what follows:
[NcA,DcA] = residue(resA,polA,[]);
polA = roots(real(DcA)), zerA = roots(real(NcA))

%%
% Now let's complicate things a little bit. A reduced order approximation, useful to simplify
% analysis and control design, is obtained by fixing a low degree, hence solving a
% least-squares problem.
% 20 Lawson iterations under the hood place our scarce resource (poles) at best, though not
% necessarily in complex conjugate pairs, so we force a recomputation. With a 2nd order
% reduction, we expect two real distinct poles:
[Hr,polAr,resAr] = aaa(GA,i*wA,'degree',2);
polAr = roots(real(poly(polAr)));
d = min(abs(i*wA(:)-polAr.'),[],1);
Q = d./(i*wA(:)-polAr.');
c = Q\GA.';
Hr = @(x) [d./(x(:)-polAr.')]*c;
[NAr] = residue(c,polAr,[]);
zerAr = roots(real(NAr)), polAr, DCgainAr = abs(Hr(0))
subplot(211), semilogx(w,20*log10(Hr(i*w)),'c-',LW,1.5), hold off
legend('G(s)','AAA','reduced order AAA',LO,SW)
subplot(212), semilogx(w,angle(Hr(i*w))*180/pi,'c-',LW,1.5), hold off
legend('G(s)','AAA','reduced order AAA',LO,SW)

%%
% To see how good AAA-LS actually is, consider the scalar example with noise found in [1],
% i.e. $f(z) = (z-1)/(z^2+z+2)$. The function is sampled at 500 logarithmically spaced points
% in the interval [0.1 10], and then normally distributed noise with a standard deviation
% of $10^{-2}$ is added:
Nc = [1 -1]; Dc = [1 1 2];
N = @(z) Nc*z.^[1 0]';
D = @(z) Dc*z.^[2 1 0]';
G = @(z) N(z(:))./D(z(:));
w = logspace(-1,1,500); magn = abs(G(i*w)); phn = -angle(G(i*w));
magn = magn+0.01*randn(1,length(magn)); phn = phn+0.01*randn(1,length(phn));
subplot(211), semilogx(w,20*log10(magn),'r-',LW,1.5), grid on, title('Magnitude (dB)')
subplot(212), semilogx(w,phn*180/pi,'r-',LW,1.5), grid on, title('Phase (degrees)')

%%
% We compute a rational approximant of degree only 2 using the above method. The AAA approximant
% shows no significant deviations from the measurements, at least in the eyeball norm. This
% time the number of Lawson iterations is
% increased, to enhance their filtering effect:
wn = [-fliplr(w) w]; magn = [fliplr(magn) magn]; phn = [-fliplr(phn) phn];
Gn = magn.*exp(i*phn);
[~,poln,resn] = aaa(Gn,i*wn,'degree',2,'lawson',30);
poln(find(real(poln)>0)) = -1;   % force system stability
poln = roots(real(poly(poln)));
dn = min(abs(i*wn(:)-poln.'),[],1);
Qn = dn./(i*wn(:)-poln.');
cn = Qn\Gn.';
Hn = @(x) [dn./(x(:)-poln.')]*cn;
 %[Nn] = residue(cn,poln,[]); zern = roots(real(Nn)), poln
 %DCgainn = abs(Hn(0))
subplot(211), hold on, semilogx(w,20*log10(Hn(i*w)),'b-',LW,1.5), hold off
legend('Noisy data','AAA approximant',LO,SW)
subplot(212), hold on, semilogx(w,angle(Hn(i*w))*180/pi,'b-',LW,1.5), hold off
legend('Noisy data','AAA approximant',LO,SW)

%%
% The poles are decently approximated, even in these perturbed conditions, as shown by the
% coefficients of the denominator:
Dcn = poly(poln)

%%
% The AAA-LS approximant effectively estimates the additive noise rather accurately,
% both on magnitude and phase:
subplot(211), loglog(wn,abs(magn(:)-abs(Hn(i*wn))),'r-',LW,1), grid on
title('Estimated noise in magnitude')
subplot(212), loglog(wn,abs(phn(:)-angle(Hn(i*wn))),'r-',LW,1), grid on
title('Estimated noise in phase')

%%
% [1] I. V. Gosea and S. G&uuml;ttel,
% Algorithms for the rational approximation of matrix-valued functions,
% arXiv:2003.06410v2, 2021.
