close all; clear; clc
% Import audio file
%file = dir('Music/Part_1/*.wav');
% [y_orig,Fs_orig] = audioread('Music/Part_1/training/DC_2.wav');
% Fs = Fs_orig/2;
% y = resample(y_orig,Fs,Fs_orig); % Subsamples signal at half the frequency

% [test_spec,tslide,n,ks] = audio_to_spec(y,Fs);
% test_spec = reshape(test_spec,n,length(tslide));
% figure()
% pcolor(tslide,ks,test_spec), 
% shading interp 
% set(gca,'Ylim',[-2000 2000]) 
% colormap(hot)
% title('Test spectrogram')

% Train classifier
[CL_data,Fs] = group_music('Music/Part_3/training/','CL');
[CL_spec,~,~,~] = audio_to_spec(CL_data,Fs);
[FK_data,Fs] = group_music('Music/Part_3/training/','FK');
[FK_spec,~,~,~] = audio_to_spec(FK_data,Fs);
[RK_data,Fs] = group_music('Music/Part_3/training/','RK');
[RK_spec,tslide,n,ks] = audio_to_spec(RK_data,Fs);
n_features = 20;
[U,S,V,CL_proj,FK_proj,RK_proj,w,cluster] = music_trainer(CL_spec,FK_spec,RK_spec,n_features);

% [DC_data,Fs] = group_music('Music/Part_1/training/','DC');
% [DC_spec,~,~,~] = audio_to_spec(DC_data,Fs);
% [DS_data,Fs] = group_music('Music/Part_1/training/','DS');
% [DS_spec,~,~,~] = audio_to_spec(DS_data,Fs);
% [KE_data,Fs] = group_music('Music/Part_1/training/','KE');
% [KE_spec,tslide,n,ks] = audio_to_spec(KE_data,Fs);
% n_features = 20;
% [U,S,V,DC_proj,DS_proj,KE_proj,w,cluster] = music_trainer(DC_spec,DS_spec,KE_spec,n_features);
%% Test classifier
[test3_CL,~] = group_music('Music/Part_3/testing/','CL');
[test3_FK,~] = group_music('Music/Part_3/testing/','FK');
[test3_RK,Fs] = group_music('Music/Part_3/testing/','RK');
test3_data = [test3_CL test3_FK test3_RK];
[test3_spec,tslide,n,ks] = audio_to_spec(test3_data,Fs);

test3_proj = U'*test3_spec;
pos = w'*test3_proj;

% classify test data
pred_class = zeros(1,length(pos));
for j = 1:length(pos)
    xy = pos(:,j);
    xy_diff = (xy-cluster(1:2,:))./cluster(3:4,:);
    [~,ind] = min(sum(abs(xy_diff),1));
    pred_class(j) = ind;
end

% Accuracy score
labels = [1 1 1 1 2 2 2 2 3 3 3 3];
correct = labels-pred_class; % 0 if correct
acc_score = 1-mean(correct~=0);
% Plot of projections on LDA basis
figure
plot(CL_proj(1,:),CL_proj(2,:),'ko')
hold on
plot(FK_proj(1,:),FK_proj(2,:),'ro')
plot(RK_proj(1,:),RK_proj(2,:),'bo')
plot(pos(1,:),pos(2,:),'k*')
xlabel('w_1')
ylabel('w_2')
legend('Classical','Folk','Rock','Test','Location','NE')
title('Test 3 data projected onto LDA basis')
%title("Projection onto LDA basis for " + n_features + " features")
hold off
% Plot energy to determine number of features to use
% sigmas = diag(S);
% 
% figure()
% subplot(1,2,1)
% plot(sigmas,'ko')
% title('Energy of singular values')
% subplot(1,2,2)
% sum_energy = cumsum(sigmas)/sum(sigmas);
% plot(sum_energy,'ko')
% title('Cumulative fraction of energy')

% Plot first 4 principal components
% for j = 1:4
%     subplot(2,2,j)
%     ut1 = reshape(U(:,j),n,length(tslide));
%     pcolor(tslide,ks,ut1), 
%     shading interp 
%     set(gca,'Ylim',[-2000 2000]) 
%     colormap(hot)
%     xlabel('Time (s)')
%     ylabel('Frequency (\omega)')
% end

% plot((1:length(v))/Fs,v);
% xlabel('Time [sec]');
% ylabel('Amplitude');
% title('Signal of Interest, v(n)');

%[ks,yt,t,tslide,ygt_spec] = audio_to_spec(y,Fs);
% FFT plot
% figure()
% subplot(2,1,1)
% plot(t(1:length(y)),y,'k')
% xlabel('Time (s)')
% ylabel('Signal')
% 
% subplot(2,1,2)
% plot(ks,abs(fftshift(yt))/max(abs(yt)),'r'); axis([-2000 2000 0 1])
% xlabel('frequency (\omega)')
% ylabel('FFT')

% Spectrogram plot
% spec1 = reshape(DC_spec(:,1),n,length(tslide));
% figure()
% pcolor(tslide,ks,spec1), 
% shading interp 
% set(gca,'Ylim',[-2000 2000]) 
% colormap(hot)
% colorbar
% title('Spectrogram of DC\_1 audio clip')
% xlabel('Time (s)')
% ylabel('Frequency (\omega)')




