function [U,S,V,v_music_1,v_music_2,v_music_3,w,clusters] = music_trainer(music_1,music_2,music_3,feature)
% Takes 3 sets of music in the form of reshaped spectrograms and
% trains a model to distinguish them.
n1 = length(music_1(1,:));
n2 = length(music_2(1,:));
n3 = length(music_3(1,:));

[U,S,V] = svd([music_1 music_2 music_3],'econ');
U = U(:,1:feature);
music = S*V'; % Projection onto principal components
music_1_proj = music(1:feature,1:n1);
music_2_proj = music(1:feature,n1+1:n1+n2);
music_3_proj = music(1:feature,n1+n2+1:n1+n2+n3);

mean_1 = mean(music_1_proj,2); % column vector with mean of each row
mean_2 = mean(music_2_proj,2);
mean_3 = mean(music_3_proj,2);

Sw = 0; % Within-class variance
for j = 1:n1
    Sw = Sw + (music_1_proj(:,j)-mean_1)*(music_1_proj(:,j)-mean_1)';
end
for j = 1:n2
    Sw = Sw + (music_2_proj(:,j)-mean_2)*(music_2_proj(:,j)-mean_2)';
end
for j = 1:n3
    Sw = Sw + (music_3_proj(:,j)-mean_3)*(music_3_proj(:,j)-mean_3)';
end

mu = (mean_1+mean_2+mean_3)/3;
Sb = n1*(mean_1-mu)*(mean_1-mu)' + n2*(mean_2-mu)*(mean_2-mu)' ...
    + n3*(mean_3-mu)*(mean_3-mu)'; % Between-class variance

% Linear discriminant analysis
[V2,D] = eig(Sb,Sw);
d = diag(D);
[~,I] = sort(d,'descend'); % Find 2 highest eigenvalues
ind1 = I(1);
ind2 = I(2);
w1 = V2(:,ind1);
w2 = V2(:,ind2);
w = [w1 w2]; % Construct plane w to project data onto
w = w/norm(w,2);

v_music_1 = w'*music_1_proj;
v_music_2 = w'*music_2_proj;
v_music_3 = w'*music_3_proj;

% Determine mean and variance of clusters
mean_v1 = mean(v_music_1,2);
mean_v2 = mean(v_music_2,2);
mean_v3 = mean(v_music_3,2);
var_v1 = var(v_music_1,0,2);
var_v2 = var(v_music_2,0,2);
var_v3 = var(v_music_3,0,2);

% Collect means and variances in matrix to output
clusters = [mean_v1 mean_v2 mean_v3;
            var_v1 var_v2 var_v3];

end

