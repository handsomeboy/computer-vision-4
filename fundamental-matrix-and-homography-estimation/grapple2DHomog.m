%% File: grappleFmatrix
%% A3 2016 handout code
%% Uses RANSAC to estimate F matrix from corresponding points.
%%
%% ADJ

clear
% close all
FALSE = 1 == 0;
TRUE = ~FALSE;
global matlabVisRoot

%% We need to ensure the path is set for the iseToolbox.
if length(matlabVisRoot)==0
  dir = pwd;
  cd C:\Users\rsanj\Dropbox\CSC2503H-Computer-Vision\Assignments\matlabVisTools\matlab
  startup;
  cd(dir);
end

reconRoot = 'C:/Users/rsanj/Dropbox/CSC2503H-Computer-Vision/Assignments/A3';
addpath([reconRoot '/data/wadham']);
addpath([reconRoot '/utils']);
addpath([reconRoot '/other_code']);


% Random number generator seed:
seed = round(sum(1000*clock));
rand('state', seed);
seed0 = seed;
% We also need to start randn. Use a seedn generated from seed:
seedn = round(rand(1,1) * 1.0e+6);
randn('state', seedn);

nTrial = 10;  % Number of ransac trials to try.

%% Wadham left image: use  wadham/001.jpg
imPath = 'data/wadham/'; fnameLeft = '001'; 
im = imread([imPath fnameLeft],'jpg');
im = double(im(:,:,2));
imDwn = blurDn(im,1)/2;
imLeft = imDwn;

%% Read correspondence data
load data/wadham/corrPnts5
%% Wadham right image data/wadham002-5.jpg use for corrPnts2-5 respectively
fnameRight = '005';
im = imread([imPath fnameRight],'jpg');
im = double(im(:,:,2));
imDwn = blurDn(im,1)/2;
imRight = imDwn;

clear imPts;
imPts = cat(3, im_pos1', im_pos2');
nPts = size(imPts,2);
if size(imPts,1) == 2
  imPts = cat(1, imPts, ones(1,nPts,2));
end

%% RANSAC for H
seeds = {};
sigma = 2.0; rho = 2;
thr=1.5;
maxInliers = 0;
% 
for kTrial = 1: nTrial
   %% Test out H matrix on a random sample of 4 points
   idTest = randperm(nPts);
   nTest = min(4, nPts);
   idTest = idTest(1:nTest);
 
   %% Solve for H matrix on the random sample
   [H HSa] = linEstH(imPts(:,idTest,1),imPts(:,idTest,2),1);
   
   %% Compute distance error
   Herror = zeros(1,nPts);
   
   % Calculate qk with Homography
   qk = H*imPts(:,:,1);
   
   % Normalize results
   qk = [qk(1,:)./qk(3,:); qk(2,:)./qk(3,:); qk(3,:)./qk(3,:)];
    
   % Compute Error Method
   Herror = imPts(:,:,2) - qk;
   Hdist = sqrt(Herror(1,:).^2 + Herror(2,:).^2);
     
   %% Detect inliers
   idInlier = Hdist < thr;
   
   %% Count inliers
   nInlier = sum(idInlier);
   if nInlier > maxInliers
     %% Store sets of sampled points with the maximum inliers
     maxInliers = nInlier;
     Hseed = struct;
     Hseed.id = idTest;
     Hseed.idInlier = idInlier;
     Hseed.nInlier = nInlier;
     Hseed.H = H;

     hSeed = length(seeds)+1;
     seeds{hSeed} = Hseed;
   end
end 

%% Best solution from RANSAC
idInlier = Hseed.idInlier;
nInlier = Hseed.nInlier;
idInlierOld = idInlier;
%% Do at most 10 iterations attempting to entrain as many points as possible.
for kIt = 1: 10
   %% Fit H using all current inliers
   [H HSa] = linEstH(imPts(:,idInlier,1), imPts(:,idInlier,2),1);
   
   %% Compute distance error
   Herror = zeros(1,nPts);
   
   %% Calculate qk with homography matrix
   qk = H*imPts(:,:,1);
   
   %% Normalize results
   qk = [qk(1,:)./qk(3,:); qk(2,:)./qk(3,:); qk(3,:)./qk(3,:)];
   
   %% Compute Error Method 2
   Herror = imPts(:,:,2) - qk;
   Hdist = sqrt(Herror(1,:).^2 + Herror(2,:).^2);
   
   %% Detect inliers
   idInlier = Hdist < thr;
   
   %% Count inliers
   nInlier = sum(idInlier);
   
   if all(idInlier == idInlierOld)
       break;
   end
   
   idInlierOld = idInlier;
end

%% Show right image
hFig = figure(5); clf;
SUPERIMPOSE = TRUE;
if SUPERIMPOSE
  image(imRight);
  colormap(gray(256));
end
resizeImageFig(hFig, size(imDwn), 1); hold on;
set(get(hFig,'CurrentAxes'),'Ydir','reverse');
hold on;
% Plot all interest point locations as blue .'s
plot(imPts(1,:,2), imPts(2,:,2), '.b');
set(gca,'YDir', 'reverse');
ax = axis;
cropBox = [ax(1) ax(3) ax(2) ax(4)];
title('Features in Right Image');
hold off

%% Show left image
hFig = figure(6); clf;
SUPERIMPOSE = true;
if SUPERIMPOSE
  image(imLeft);
  colormap(gray(256));
end
resizeImageFig(hFig, size(imDwn), 1); hold on;
set(get(hFig,'CurrentAxes'),'Ydir','reverse');
hold on;
% Plot all interest point locations as blue .'s
plot(imPts(1,:,1), imPts(2,:,1), '.b');
set(gca,'YDir', 'reverse');
ax = axis;
cropBox = [ax(1) ax(3) ax(2) ax(4)];
title('Features in Left Image');
hold off

%% Calculate warped image
warped_left_image = homogWarp(imLeft,H^-1);

%% Show warped left image
hFig = figure(7); clf;
SUPERIMPOSE = true;
if SUPERIMPOSE
  image(warped_left_image);
  colormap(gray(256));
end
resizeImageFig(hFig, size(imDwn), 1); hold on;
set(get(hFig,'CurrentAxes'),'Ydir','reverse');
hold on;
set(gca,'YDir', 'reverse');
ax = axis;
cropBox = [ax(1) ax(3) ax(2) ax(4)];
title('Warped Left Image');
hold off

%% Calculate warped image
warped_right_image = homogWarp(imRight,H);

%% Show warped left image
hFig = figure(8); clf;
SUPERIMPOSE = true;
if SUPERIMPOSE
  image(warped_right_image);
  colormap(gray(256));
end
resizeImageFig(hFig, size(imDwn), 1); hold on;
set(get(hFig,'CurrentAxes'),'Ydir','reverse');
hold on;
set(gca,'YDir', 'reverse');
ax = axis;
cropBox = [ax(1) ax(3) ax(2) ax(4)];
title('Warped Right Image');
hold off

%% Show warped left image will all points
hFig = figure(9); clf;
SUPERIMPOSE = true;
if SUPERIMPOSE
  image(warped_left_image);
  colormap(gray(256));
end
resizeImageFig(hFig, size(imDwn), 1); hold on;
set(get(hFig,'CurrentAxes'),'Ydir','reverse');
hold on;
% Plot all interest point locations as blue .'s
plot(qk(1,:),qk(2,:),'.r');
plot(imPts(1,:,1), imPts(2,:,1), '.b');

% plot([qk(1,1),imPts(1,1,1)],[qk(2,1),imPts(2,1,1)]);
plot([qk(1,:),imPts(1,:,1)],[qk(2,:),imPts(2,:,1)]);

set(gca,'YDir', 'reverse');
ax = axis;
cropBox = [ax(1) ax(3) ax(2) ax(4)];
title('Warped Left Image w/ All Mapping');
hold off

%% Show warped left image with only inliers
hFig = figure(10); clf;
SUPERIMPOSE = true;
if SUPERIMPOSE
  image(warped_left_image);
  colormap(gray(256));
end
resizeImageFig(hFig, size(imDwn), 1); hold on;
set(get(hFig,'CurrentAxes'),'Ydir','reverse');
hold on;
% Plot all interest point locations as blue .'s
plot(qk(1,idInlier),qk(2,idInlier),'.r');
plot(imPts(1,idInlier,1), imPts(2,idInlier,1), '.b');

% plot([qk(1,1),imPts(1,1,1)],[qk(2,1),imPts(2,1,1)]);
plot([qk(1,idInlier),imPts(1,idInlier,1)],[qk(2,idInlier),imPts(2,idInlier,1)]);

set(gca,'YDir', 'reverse');
ax = axis;
cropBox = [ax(1) ax(3) ax(2) ax(4)];
title('Warped Left Image w/ All Mapping');
hold off

