function [f_max,tpulse,yft] = get_frequencies(y,Fs,a,minpeakdist,minpeakheight,flag)
% Takes a vector y of signal data and sampling frequency Fs and
% returns the center frequencies in the signal

% Define domains
n = pow2(nextpow2(length(y))); % 2^20
L = n/Fs;
t2 = [0:1/Fs:L];
t = t2(1:n);
k = (2*pi/L)*[0:n/2-1 -n/2:-1];
ks = fftshift(k);

y_pad = [y zeros(1,n-length(y))]; % zero-padding y

% Find location of pulses
[pks,tpulse] = findpeaks(y_pad,Fs,'MinPeakDistance',minpeakdist,'MinPeakHeight',minpeakheight);
yft = zeros(length(tpulse),n);
plot_flag = flag; % Make plots if true
for ii = 1:length(tpulse)
    g = exp(-a*(t-tpulse(ii)).^2);
    yf = g.*y_pad;
    yft(ii,:) = fft(yf);
    if plot_flag
        figure()
        subplot(3,1,1), plot(t,y_pad,'k',t,g,'r')
        subplot(3,1,2), plot(t,yf,'k')
        subplot(3,1,3), plot(ks,abs(fftshift(yft(ii,:)))/max(abs(yft(ii,:))))
        axis([-20000 20000 0 1])
    end
end

% Find maxima of frequency data
[M,I] = max(abs(yft),[],2);
omega = k(I);
f_max = omega/(2*pi);

end

