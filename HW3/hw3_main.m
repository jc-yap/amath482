close all; clear; clc
myFiles = dir('vids\*.mat');
test1 = zeros(6,181);
test2 = zeros(6,241);
test3 = zeros(6,151);
test4 = zeros(6,281);
count = 1;
for j = 1:length(myFiles)
    file = strcat('vids\',myFiles(j).name);
    pos_vec = mass_position(file);
%     figure()
%     subplot(3,1,1)
%     plot(pos_vec(:,1),'k.')
%     subplot(3,1,2)
%     plot(pos_vec(:,2),'k')
%     subplot(3,1,3)
%     plot(pos_vec(:,1),pos_vec(:,2),'k')
    
    if mod(j,4) == 1
        test1(count:count+1,:) = pos_vec';
    elseif mod(j,4) == 2
        test2(count:count+1,:) = pos_vec';
    elseif mod(j,4) == 3
        test3(count:count+1,:) = pos_vec';
    else
        test4(count:count+1,:) = pos_vec';
    end
        
    if mod(j,4) == 0
        count = count+2;
    end
end

save('test1.mat','test1')
save('test2.mat','test2')
save('test3.mat','test3')
save('test4.mat','test4')
%%
close all; clear all; clc
% subtract mean
load('test1.mat')
load('test2.mat')
load('test3.mat')
load('test4.mat')
test = test4;
[m,n] = size(test);
mn = mean(test,2);
X1 = test - repmat(mn,1,n);

[U,S,V] = svd(X1,'econ');
% Compute and plot rank-1 approximations
X_rank1 = U(:,1)*S(1,1)*V(:,1).';
figure()
plot(X_rank1(1,:),X_rank1(2,:),'k.');
axis equal
hold on
plot(X_rank1(3,:),X_rank1(4,:),'r.');
plot(X_rank1(5,:),X_rank1(6,:),'b.');
legend('Camera 1','Camera 2','Camera 3')
hold off
% saveas(gcf,'test4_1pc.png')
% Plot first principal component
y1 = (S(1,1)/sqrt(n-1))*U(:,1);
figure()
compass(y1(1),y1(2));
hold on
axis equal

% Plot second principal component
y2 = (S(2,2)/sqrt(n-1))*U(:,2);
compass(y2(1),y2(2));
hold off
Y = U'*test2;
%plot(Y(1,:),Y(2,:),'k.')

% Plot energy
figure()
subplot(1,2,1)
sig = diag(S);
plot(sig,'ko')
title('Energy of singular values')
subplot(1,2,2)
sum_energy = cumsum(sig)/sum(sig);
plot(sum_energy,'ko')
title('Cumulative fraction of energy')
% saveas(gcf,'test4_energy.png')