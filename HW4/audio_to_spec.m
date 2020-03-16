function [music_spec,tslide,n,ks] = audio_to_spec(y_mat,Fs)
% Takes in an audio signal and returns a matrix to plot a spectrogram
[rows,cols] = size(y_mat);
n = pow2(nextpow2(rows)); % number of fourier modes
L = n/Fs;

t2 = linspace(0,L,n+1); % time discretization
t = t2(1:n);

k = (2*pi/L)*[0:n/2-1 -n/2:-1]; % frequency components
ks = fftshift(k);
%yt = fft(y,n); % Fourier transforms y with zero-padding

% Construct Gabor window
tslide = 0:0.1:5; % where the window is centered
a = 100; % controls window width

% Initialize return matrix
music_spec = zeros(length(tslide)*n,cols);
for ii = 1:cols
    ygt_spec = zeros(n,length(tslide));
    y = y_mat(:,ii);
    for j=1:length(tslide)
        g=exp(-a*(t(1:rows)-tslide(j)).^2);
        yg=g'.*y; 
        ygt=fft(yg,n);
        ygt_spec(:,j) = rescale(fftshift(abs(ygt))); % Scaled
    end
    music_spec(:,ii) = reshape(ygt_spec,length(tslide)*n,1);

end

