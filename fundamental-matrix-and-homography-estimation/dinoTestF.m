%% File: dinoTestF
%% A3 2016 handout code
%% Estimate F matrix from synthetic corresponding points.
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
addpath([reconRoot '/utils']);

% Random number generator seed:
seed = round(sum(1000*clock));
rand('state', seed);
seed0 = seed;
% We also need to start randn. Use a seedn generated from seed:
seedn = round(rand(1,1) * 1.0e+6);
randn('state', seedn);

%% Loop parameters
nTrial = 10;  %% Number of ransac trials to use
nNoiseTrial = 25; %% Number of noise trials to use
nSigmaTrial = 5/0.1; %% Number of sigma values to use

%% Global Scaling Parameters
NUM_RESCALE = 1; %% 1 - Default, 0 - Check
SCLZ = 1; %% 1 - Default, 0 - Planar, 0.1 - Flat (Check)
D_SCALE_FACTOR = 1; %% 1 - Default, 0.1 - Check

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up cameras
%% The cameras are automatically rotated by projectDino to fixate
%% on the mean of the 3D points.  We do not always want to allow
%% the cameras to move in this fashion (eg. when we change sclZ).
%% So we will compute the rotations for the left and right cameras
%% once and for all, and then use these.

f = 100; % focal length
dLeft = [-50*D_SCALE_FACTOR, 0, -150]';  % Center of projection for left camera
dRight = [50*D_SCALE_FACTOR, 0, -150]';  % Center of projection for right camera
%% Compute camera rotations to fixate on Dino's center.
[pLeft polys MintLeft MextLeft] = projectDino(f, dLeft, [], 1.0);
Rleft = MextLeft(:, 1:3);
[pRight polys MintRight MextRight] = projectDino(f, dRight, [], 1.0);
Rright = MextRight(:, 1:3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Generate data...
sclZ = SCLZ;
%% Dino left image data
[pLeft polys MintLeft MextLeft] = projectDino(f, dLeft, Rleft, sclZ);

%% Dino right image data
[pRight polys MintRight MextRight] = projectDino(f, dRight, Rright, sclZ);

%% Show left and right images
hFig = figure(1); clf; 
plot(pLeft(1,:), pLeft(2,:), '.b');
axis xy; axis equal;
xlabel('X'); ylabel('Y');
axis([-150 150 -100 100]);
title('Left image of Dino');
pause(0.1);

hFig = figure(2); clf; 
plot(pRight(1,:), pRight(2,:), '.b');
axis xy; axis equal;
axis([-150 150 -100 100]);
xlabel('X'); ylabel('Y');
title('Right image of Dino');
pause(0.1);

%% Compute Ground Truth
[F0] = findFGroundTruth(MintLeft,MextLeft,MintRight,MextRight);

% Sanity check
check = pLeft'*F0*pRight;
check = diag(check);
F0check = max(max(check))

%% Build correspondence data
clear imPts;
imPts = cat(3, pLeft, pRight);
nPts = size(imPts,2);
if size(imPts,1) == 2
  imPts = cat(1, imPts, ones(1,nPts,2));
end

%% RANSAC for F
seeds = {};
sigma = 2.0; rho = 2;
for kTrial = 1: nTrial
  %% Test out F matrix on a random sample of 8 points
  idTest = randperm(nPts);
  nTest = min(8, nPts);
  idTest = idTest(1:nTest);

  %% Solve for F matrix on the random sample
  [F Sa Sf] = linEstF(imPts(:,idTest,1), imPts(:,idTest,2),1);

  %% Compute perpendicular error of all points to epipolar lines
  perpErrL = zeros(1,nPts);
  for k = 1:nPts
    lk = imPts(:,k,2)' * F';
    perpErrL(k) = (lk* imPts(:,k,1))/norm(lk(1:2));
  end

  %% Detect inliers
  idInlier = abs(perpErrL) < rho*sigma;

  %% Count inliers
  nInlier = sum(idInlier);
  if nInlier > 20
    %% Store sets of sampled points with at least 20 inliers
    seed = struct;
    seed.id = idTest;
    seed.idInlier = idInlier;
    seed.nInlier = nInlier;
    seed.F = F;

    kSeed = length(seeds)+1
    seeds{kSeed} = seed;
  end
end 
%% Done RANSAC trials

%% Extract best solution
nInliers = zeros(1, length(seeds));
for ks = 1:length(seeds)
  nInliers(ks) = seeds{ks}.nInlier;
end 
[nM ks] = max(nInliers);
nInliers(ks)

%%  Refine estimate of F using all inliers.
F = seeds{ks}.F;
idInlier = seeds{ks}.idInlier;

idInlierOld = idInlier;
sum(idInlier)
%% Do at most 10 iterations attempting to entrain as many points as possible.
for kIt = 1: 10
  %% Fit F using all current inliers
  [F Sa Sf] = linEstF(imPts(:,idInlier,1), imPts(:,idInlier,2),1);

  %% Compute perpendicular error to epipolar lines
  perpErrL = zeros(1,nPts);
  for k = 1:nPts
    lk = imPts(:,k,2)' * F';
    perpErrL(k) = (lk* imPts(:,k,1))/norm(lk(1:2));
  end
  idInlier = abs(perpErrL) < rho*sigma;
  nInlier = sum(idInlier)

  %% If we have the same set of inliers as the previous iteration then stop.
  if all(idInlier == idInlierOld)
    break;
  end
  idInlierOld = idInlier;
end

%%%%%%%%%%%%%%%%%%%%% Plot results
nTest = 64;  %% Number of epipolar lines to plot
nCol = 16;   %% Number of different colours to use.
col = hsv(nCol);  %% Colour map.

%% Random sample the lines to plot
idLines = randperm(nPts);  
idLines = idLines(1:nTest);

%% Show left image
hFig = figure(3);
clf; hold on;
% Plot all interest point locations as blue .'s
plot(imPts(1,:,1), imPts(2,:,1), '.b');
axis([-150 150 -100 100]); axis xy; axis equal;
ax = axis;
cropBox = [ax(1) ax(3) ax(2) ax(4)];
title('Epipolar Lines in Left Image');
for kl = 1:length(idLines)
  % Plot interest point location corresponding to epipolar line as a "o"
  % in the same colour as the epipolar line.
  k = idLines(kl);
  plot(imPts(1,k,1), imPts(2,k,1), 'o', 'Color', col(mod(k,nCol)+1,:));
  % Plot epipolar line.
  lk = imPts(:,k,2)' * F';
  epk = cropLineInBox(lk(1:2), lk(3), cropBox); 
  set(line(epk(:,1), epk(:,2)), 'Color', col(mod(k,nCol)+1,:));
end

%% Show right image
hFig = figure(4);
clf; hold on;
% Plot all interest point locations as blue .'s
plot(imPts(1,:,2), imPts(2,:,2), '.b');
axis([-150 150 -100 100]); axis xy; axis equal;
ax = axis;
cropBox = [ax(1) ax(3) ax(2) ax(4)];
title('Epipolar Lines in Right Image');
perpErrR = [];
for kl = 1:length(idLines)
  % Plot interest point location corresponding to epipolar line as a "o"
  % in the same colour as the epipolar line.
  k = idLines(kl);
  plot(imPts(1,k,2), imPts(2,k,2), 'o', 'Color', col(mod(k,nCol)+1,:));
  % Plot epipolar line.
  lk = imPts(:,k,1)' * F;
  epk = cropLineInBox(lk(1:2), lk(3), cropBox); 
  set(line(epk(:,1), epk(:,2)), 'Color', col(mod(k,nCol)+1,:));
  perpErrR = [perpErrR (lk * imPts(:,k,2)/norm(lk(1:2)))];
end

%% Compute perpendicular distance to epipolar lines in left and right images.
perpErrL = [];
for k = 1:nPts
  lk = imPts(:,k,2)' * F';
  perpErrL = [perpErrL (lk * imPts(:,k,1))/norm(lk(1:2))];
end
perpErrR = [];
for k = 1:nPts
  lk = imPts(:,k,1)' * F;
  perpErrR = [perpErrR (lk * imPts(:,k,2)/norm(lk(1:2)))];
end

%% Plot a histogram of the perpendicular distances
err = [perpErrL perpErrR];
err = min(err, 10);
err = max(err, -10);
[n b] = histo(err, 64);
figure(5); clf;
plot(b,n);
title('Distance to epipolar line');
xlabel('Error in pixels');
ylabel('Frequency');

%% Count inliers
idL = abs(perpErrL)< rho*sigma;
idR = abs(perpErrR) < rho*sigma;
idInlier = idL & idR;
sum(idInlier)
sum(idInlier)/nPts

%% save 'Fcorr' F Sa Sf idInlier nInliers

%% Sanity check for F
check2 = pLeft'*F*pRight;
check2 = diag(check2);
Fcheck = max(max(check2))

%% Part 2 - Find Perpendicular Error
region = [-70 0 -40 30]; %% [X0 Y0 X1 Y1]
[FError] = findErrorInFEst(pRight, region, F,F0)

%% Part 3 - Add Noise
sigmaN = 0;
medianFError = zeros(nSigmaTrial,1);
sigmaFError = zeros(nSigmaTrial,1);

pLeftwNoise(3,:) = pLeft(3,:);
pRightwNoise(3,:) = pRight(3,:);

for kSigmaTrial = 1:nSigmaTrial
    sigmaFError(kSigmaTrial) = sigmaN;
    Ferrors = zeros(nNoiseTrial,1);

    for kNoiseTrial = 1:nNoiseTrial
        % Add noise     
        pLeftwNoise(1:2,:) = pLeft(1:2,:) + sigmaN.*randn(2,nPts);
        pRightwNoise(1:2,:) = pRight(1:2,:) + sigmaN.*randn(2,nPts);

        % Find F
        [F2 Sa Sf] = linEstF(pLeftwNoise, pRightwNoise, NUM_RESCALE);

        % Get Error
        [Ferrors(kNoiseTrial)] = findErrorInFEst(pRight,region, F2,F0);

    end

    % Get Median Error
    medianFError(kSigmaTrial) = median(Ferrors);

    % Increment Sigma
    sigmaN = sigmaN + 0.1;

end

%% Plot sigma vs. errors

hFig = figure(6);clf;
hold on
if NUM_RESCALE == 0 && SCLZ == 1 && D_SCALE_FACTOR == 1
    title('SigmaN vs. Median Error w/ Hartley Rescaling = 0');
elseif NUM_RESCALE == 1 && SCLZ == 0.1 && D_SCALE_FACTOR == 1
    title('SigmaN vs. Median Error w/ sclZ = 0.1');
elseif NUM_RESCALE == 1 && SCLZ == 1 && D_SCALE_FACTOR == 0.1
    title('SigmaN vs. Median Error w/ Camera Distance = 10');
elseif NUM_RESCALE == 1 && SCLZ == 1 && D_SCALE_FACTOR == 1
    title('SigmaN vs. Median Error');
else
    title('SigmaN vs. Median Error w/ invalid case');
end
scatter(sigmaFError,medianFError,'filled');
hold off

hFig = figure(7);clf;
hold on
if NUM_RESCALE == 0 && SCLZ == 1 && D_SCALE_FACTOR == 1
    title('SigmaN vs. Median Error w/ Hartley Rescaling = 0');
elseif NUM_RESCALE == 1 && SCLZ == 0.1 && D_SCALE_FACTOR == 1
    title('SigmaN vs. Median Error w/ sclZ = 0.1');
elseif NUM_RESCALE == 1 && SCLZ == 1 && D_SCALE_FACTOR == 0.1
    title('SigmaN vs. Median Error w/ Camera Distance = 10');
elseif NUM_RESCALE == 1 && SCLZ == 1 && D_SCALE_FACTOR == 1
    title('SigmaN vs. Median Error');
else
    title('SigmaN vs. Median Error w/ invalid case');
end
plot(sigmaFError,medianFError)
xlabel('\sigma_n');
ylabel('Median Error');
hold off





