function [result_mat,Fs] = group_music(sub_dir,artist)
% Looks in the specified directory and groups .wav files of the
% specified type into matrices for training and testing. Subsamples
% the given audio at half the frequency.
% n_files = length(dir(sub_dir));
if contains(sub_dir,'training')
    range_file = [1:16];
else
    range_file = [17:20];
end
result_mat = zeros(110250,length(range_file));
for j = range_file
    file_name = strcat(sub_dir,artist,'_',num2str(j),'.wav');
    [y_orig,Fs_orig] = audioread(file_name);
    Fs = Fs_orig/2;
    y = resample(y_orig,Fs,Fs_orig); % Subsampling
    result_mat(:,mod(j-1,16)+1) = y;
end    

end

