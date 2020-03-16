function [pos_vec] = mass_position(filename)
% Takes in a .mat video file and returns a matrix of the
% positions of the mass in each frame. Column 1 returns
% x-position and column 2 returns y_position.
varName = strcat('vidFrames',filename(9:11));
data = load(filename,varName);
vidFrames = data.(varName);

% load('cam2_4.mat');
% load('cam3_4.mat');
numFrames = size(vidFrames,4);
vid_gs = zeros(size(vidFrames,[1 2 4]));
% Convert to grayscale
for j = 1:numFrames
    X = vidFrames(:,:,:,j);
    X = rgb2gray(X);
    vid_gs(:,:,j) = X;
end
%imtool(vid_gs(:,:,1),[])
dl = 46;
if strcmp('vidFrames1_1',varName)
    rstart = 266;
    cstart = 323;
    frame_start = 30;
    n_frames = 180;
elseif strcmp('vidFrames1_2',varName)
    rstart = 317;
    cstart = 317;
    frame_start = 13;
    n_frames = 240;
elseif strcmp('vidFrames1_3',varName)
    rstart = 309;
    cstart = 321;
    frame_start = 59;
    n_frames = 150;
elseif strcmp('vidFrames1_4',varName)
    rstart = 287;
    cstart = 373;
    frame_start = 72;
    n_frames = 280;
elseif strcmp('vidFrames2_1',varName)
    rstart = 297;
    cstart = 267;
    frame_start = 40;
    n_frames = 180;
elseif strcmp('vidFrames2_2',varName)
    rstart = 329;
    cstart = 286;
    dl = 2*dl;
    frame_start = 38;
    n_frames = 240;
elseif strcmp('vidFrames2_3',varName)
    rstart = 300;
    cstart = 221;
    dl = 2*dl;
    frame_start = 86;
    n_frames = 150;
elseif strcmp('vidFrames2_4',varName)
    rstart = 262;
    cstart = 231;
    dl = 1.5*dl;
    frame_start = 79;
    n_frames = 280;
elseif strcmp('vidFrames3_1',varName)
    rstart = 268;
    cstart = 332;
    frame_start = 29;
    n_frames = 180;
elseif strcmp('vidFrames3_2',varName)
    rstart = 249;
    cstart = 350;
    frame_start = 18;
    n_frames = 240;
    % cut out last few frames
elseif strcmp('vidFrames3_3',varName)
    rstart = 224;
    cstart = 365;
    frame_start = 50;
    n_frames = 150;
    % cut out last few frames
elseif strcmp('vidFrames3_4',varName)
    rstart = 180;
    cstart = 345;
    dl = 2*dl;
    frame_start = 71;
    n_frames = 280;
end
can = vid_gs(rstart:rstart+dl,cstart:cstart+dl);
pos_vec = zeros(numFrames,2);

for j = 1:numFrames
    framej = vid_gs(:,:,j);
    c = normxcorr2(can,framej);
    % figure, surf(c), shading flat
    
    if j == 1
        [ypeak, xpeak] = find(c==max(c(:)));
    else
        c_red = c(ypeak-dl:ypeak+dl,xpeak-dl:xpeak+dl);
        max_red = max(c_red(:));
        [ypeak, xpeak] = find(c==max_red);
    end
    yoffSet = ypeak-size(can,1);
    xoffSet = xpeak-size(can,2);
%     imshow(framej,[]);
%     hold on
%     drawrectangle('Position',[xoffSet+1, yoffSet+1, dl+1, dl+1]);
%     scatter(xoffSet+1+(dl+1)/2,yoffSet+1+(dl+1)/2,'y.')
%     drawnow
%     hold off
    pos_vec(j,:) = [xoffSet+1+(dl+1)/2,yoffSet+1+(dl+1)/2];
    can = framej(yoffSet+1:yoffSet+dl+2,xoffSet+1:xoffSet+dl+2);
end

pos_vec = pos_vec(frame_start:frame_start+n_frames,:);
end

