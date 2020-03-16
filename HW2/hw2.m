close all; clear; clc
load handel
v = y';
t = (1:length(v))/Fs;
subplot(2,1,1)
plot(t,v)
xlim([0 8.5])
hold on
xlabel('Time [sec]');
ylabel('Amplitude');
title('Signal of Interest, v(n)');

% Play audio
% p8 = audioplayer(v,Fs);
% playblocking(p8);

% Define domains
n = pow2(nextpow2(length(v))); % 2^17
L = n/Fs;
t2 = [0:1/Fs:L];
t = t2(1:n);
k = (2*pi/L)*[0:n/2-1 -n/2:-1];
ks = fftshift(k);

% Define Gabor filter
a = 2; % controls filter width
tau = 2; % center of filter
filter = exp(-a*(t-tau).^2);
plot(t,filter,'Linewidth',2)
hold off
legend('Signal','Filter')

% Apply Gabor filter to signal
v_pad = [v zeros(1,n-length(v))]; % zero-padding v
v_filt = filter.*v_pad;
v_filt_t = fft(v_filt);

subplot(2,1,2)
plot(t,v_filt)
xlim([0 8.5])
xlabel('Time [sec]');
ylabel('Amplitude');
title('Filtered signal');
saveas(gcf,'gab_filter.png')
% figure()
% plot(ks,abs(fftshift(v_filt_t)))

% figure()
% Slide Gabor window across time domain
tslide = [0:0.1:t(length(v))];
vft_spec = zeros(length(tslide),n);
plot_flag = false; % Plot animation if true
for ii = 1:length(tslide)
    g = exp(-a*(t-tslide(ii)).^2);
    vf = g.*v_pad;
    vft = fft(vf);
    vft_spec(ii,:) = fftshift(abs(vft));
    if plot_flag == true
        subplot(3,1,1), plot(t,v_pad,'k',t,g,'r')
        subplot(3,1,2), plot(t,vf,'k')
        subplot(3,1,3), plot(ks,abs(fftshift(vft))/max(abs(vft)))
        axis([-20 20 0 1])
        drawnow
        pause(0.1)
    end
end

% Plot spectrogram
% pcolor(tslide,ks,vft_spec.'), shading interp
% set(gca, 'Ylim',[-20,20],'Fontsize',14)
% colormap(hot)
% colorbar
% caxis([0 200])

% Spectrograms for varying window sizes
figure()
a_vec = [10 5 2 0.1];
for jj = 1:length(a_vec)
    a = a_vec(jj);
    tslide = [0:0.1:t(length(v))];
    vft_spec = zeros(length(tslide),n);
    for ii = 1:length(tslide)
        g = exp(-a*(t-tslide(ii)).^2);
        vf = g.*v_pad;
        vft = fft(vf);
        vft_spec(ii,:) = fftshift(abs(vft));
    end
    subplot(2,2,jj)
    pcolor(tslide,ks,vft_spec.'), shading interp
    set(gca, 'Ylim',[-20,20],'Fontsize',12)
    title(['a = ',num2str(a)],'Fontsize',12)
    xlabel('Time (s)')
    ylabel('Frequency')
    colormap(hot)
    caxis([0 100])
end
saveas(gcf,'spec_1.png')

%% Explore over and undersampling
dt_vec = [0.1 0.5 1 5];
for jj = 1:length(dt_vec)
    a = 2;
    dt = dt_vec(jj);
    tslide = [0:dt:t(length(v))];
    vft_spec = zeros(length(tslide),n);
    for ii = 1:length(tslide)
        g = exp(-a*(t-tslide(ii)).^2);
        vf = g.*v_pad;
        vft = fft(vf);
        vft_spec(ii,:) = fftshift(abs(vft));
    end
    subplot(2,2,jj)
    pcolor(tslide,ks,vft_spec.'), shading interp
    set(gca, 'Ylim',[-20,20],'Fontsize',12)
    title(['dt = ',num2str(dt)],'Fontsize',12)
    xlabel('Time (s)')
    ylabel('Frequency')
    colormap(hot)
    caxis([0 100])
end
saveas(gcf,'spec_sampling.png')

%% Part 2
% Piano
[y,Fs] = audioread('music1.wav');
tr_piano=length(y)/Fs; % record time in seconds
% plot((1:length(y))/Fs,y);
% xlabel('Time [sec]'); ylabel('Amplitude');
% title('Mary had a little lamb (piano)');
% p8 = audioplayer(y,Fs); playblocking(p8);

a = 60;
mindist = 0.2;
minheight = 0.2;
[f_piano,t_piano,yft_piano] = get_frequencies(y',Fs,a,mindist,minheight,false);

figure()
subplot(3,1,1)
plot(ks,abs(fftshift(yft_piano(1,:)))/max(abs(yft_piano(1,:))))
axis([-10000 10000 0 1])
title('Frequencies for note E')
xlabel('Frequency')
ylabel('Magnitude')
subplot(3,1,2)
plot(ks,abs(fftshift(yft_piano(2,:)))/max(abs(yft_piano(2,:))))
axis([-10000 10000 0 1])
title('Frequencies for note D')
xlabel('Frequency')
ylabel('Magnitude')
subplot(3,1,3)
plot(ks,abs(fftshift(yft_piano(3,:)))/max(abs(yft_piano(3,:))))
axis([-10000 10000 0 1])
title('Frequencies for note C')
xlabel('Frequency')
ylabel('Magnitude')

saveas(gcf,'overtones_piano.png')
%%
% Define domains
n = pow2(nextpow2(length(y))); % 2^20
L = n/Fs;
t2 = [0:1/Fs:L];
t = t2(1:n);
k = (2*pi/L)*[0:n/2-1 -n/2:-1];
ks = fftshift(k);

figure()
% Define Gabor filter
tau_vec = [0.8 1.25]; % center of filter
for ii = 1:2
    tau = tau_vec(ii);
    filter = exp(-a*(t-tau).^2);

    % Apply Gabor filter to signal
    y_pad = [y' zeros(1,n-length(y))]; % zero-padding y
    y_filt = filter.*y_pad;
    y_filt_t = fft(y_filt);

    % Make plots
    subplot(2,1,ii)
    plot(t,y_pad)
    xlabel('Time [sec]'); ylabel('Amplitude');
    title(['Filter applied to note ',num2str(ii)]);
    axis([0 15 -1 1])
    hold on
    plot(t,filter);
    hold off
end
saveas(gcf,'piano_filter.png')

%%
% Recorder
figure(2)
[y,Fs] = audioread('music2.wav');
tr_rec=length(y)/Fs; % record time in seconds

figure()
plot((1:length(y))/Fs,y);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Mary had a little lamb (recorder)');
% p8 = audioplayer(y,Fs); playblocking(p8);
[f_rec,t_rec,yft_rec] = get_frequencies(y',Fs,60,0.3,0.08,false);

figure()
subplot(3,1,1)
plot(ks,abs(fftshift(yft_rec(1,:)))/max(abs(yft_rec(1,:))))
axis([-30000 30000 0 1])
title('Frequencies for note B')
xlabel('Frequency')
ylabel('Magnitude')
subplot(3,1,2)
plot(ks,abs(fftshift(yft_rec(2,:)))/max(abs(yft_rec(2,:))))
axis([-30000 30000 0 1])
title('Frequencies for note A')
xlabel('Frequency')
ylabel('Magnitude')
subplot(3,1,3)
plot(ks,abs(fftshift(yft_rec(3,:)))/max(abs(yft_rec(3,:))))
axis([-30000 30000 0 1])
title('Frequencies for note G')
xlabel('Frequency')
ylabel('Magnitude')

saveas(gcf,'overtones_rec.png')
%% Generate music score
notes = [3 2 1 2 3 3 3 2 2 2 3 3 3 3 2 1 2 3 3 3 3 2 2 3 2 1];
t_piano = t_piano(1:26);
t_rec = t_rec(1:26);
figure()
plot(t_piano,2+notes,'ko','MarkerFaceColor', 'k')
title('Music score of piano')
xlabel('Time (s)')
ylabel('Scale')
scale = {'A';'B';'C';'D';'E';'F';'G'};
set(gca,'ytick',[1:7],'yticklabel',scale)
grid on
ylim([1,7])
saveas(gcf,'score_piano.png')

figure()
plot(t_rec,notes,'ro','MarkerFaceColor', 'r')
title('Music score of recorder')
xlabel('Time (s)')
ylabel('Scale')
scale = {'G';'A';'B';'C';'D';'E';'F'};
set(gca,'ytick',[1:7],'yticklabel',scale)
grid on
ylim([1,7])
saveas(gcf,'score_rec.png')