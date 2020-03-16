clear; close all; clc; 
load Testdata

L=15; % spatial domain 
n=64; % Fourier modes 
x2=linspace(-L,L,n+1); 
x=x2(1:n); y=x; z=x; 
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; ks=fftshift(k);

[X,Y,Z]=meshgrid(x,y,z); 
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);

Unt_avg = zeros(1,n^3);
for j=1:20 
    Un = Undata(j,:);
    Unt = fft(Un);
    Unt_avg = Unt_avg + Unt;
end
Unt_avg = abs(fftshift(Unt_avg))./20;
[M,I] = max(Unt_avg);

% Convert indices to frequencies
Kxmax = Kx(I);
Kymax = Ky(I);
Kzmax = Kz(I);

% [Kxmax,Kymax,Kzmax] = [-4.8171,-6.7021,-1.0472]

% Construct Gaussian filter centered at [Kxmax,Kymax,Kzmax]
tau = 0.4;
filter = exp(-tau*((Kx-Kxmax).^2+(Ky-Kymax).^2+(Kz-Kzmax).^2));
filter = reshape(filter,[1,n^3]);

% Initialize position variable, stores (x,y,z) coordinates of
% marble at each time point
pos = zeros(20,3);

% Loop through time points of given data
for j = 1:20
    % Filter frequencies from each time point in given data and
    % transform back to spatial domain to plot
    Un1 = Undata(j,:);
    Un1t = fft(Un1);
    Un1t = fftshift(Un1t);
    Un1tf = filter.*Un1t;
    Un1f = ifft(ifftshift(Un1tf));
    
    % Determine position of marble at each time
    [m,ii] = max(Un1f);
    xj = X(ii);
    yj = Y(ii);
    zj = Z(ii);
    pos(j,:) = [xj yj zj];
    
    isosurface(X,Y,Z,abs(reshape(Un1f,[n,n,n])),0.2)
    hold on
    view([-1,-1,0.5])
    title('Position of marble at each time')
    xlabel('x')
    ylabel('y')
    zlabel('z')
    axis([-20 20 -20 20 -20 20]) 
    grid on
    pause(0.2)
end
print -depsc m_pos.eps
hold off

% Plot path of marble
figure()
plot3(pos(:,1),pos(:,2),pos(:,3),'LineWidth',2)
title('Trajectory of marble')
xlabel('x')
ylabel('y')
zlabel('z')
axis([-20 20 -20 20 -20 20])
grid on
print -depsc m_traj.eps

% Determine final position of marble
final_pos = pos(end,:);