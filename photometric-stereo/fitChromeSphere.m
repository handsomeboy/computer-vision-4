function [L] = fitChromeSphere(chromeDir, nDir, chatty)
  % [L] = fitChromeSphere(chromeDir, nDir, chatty)
  % Input:
  %  chromeDir (string) -- directory containing chrome images.
  %  nDir -- number of different light source images.
  %  chatty -- true to show results images. 
  % Return:
  %  L is a 3 x nDir image of light source directions.

  % Since we are looking down the z-axis, the direction
  % of the light source from the surface should have
  % a negative z-component, i.e., the light sources
  % are behind the camera.
    
  if ~exist('chatty', 'var')
    chatty = false;
  end
    
  mask = ppmRead([chromeDir, 'chrome.mask.ppm']);
  mask = mask(:,:,1) / 255.0;
  
  % Get center of circle
  %https://www.mathworks.com/help/images/examples/detect-and-measure-circular-objects-in-an-image.html
  %figure;
  %showIm(mask);
  %d = imdistline; %Used to find radius range (100-125)
  [center, r] = imfindcircles(mask,[100 125]);
  xc = center(1);
  yc = center(2);
  
  camera_dir = [0; 0; -1];
  
  for n=1:nDir
    fname = [chromeDir,'chrome.',num2str(n-1),'.ppm'];
    im = ppmRead(fname);
    imData(:,:,n) = im(:,:,1);           % red channel
    
    % mask image
    imData(:,:,n) = imData(:,:,n).*mask(:,:,1); 
    
    % Find pixel coordinate of the maximum value
    [M1,I1] = max(imData(:,:,n));
    [M2,I2] = max(M1);
    y(n) = I1(I2);
    x(n) = I2;
    %fprintf('max: (%i,%i) \n', pt_x(n), pt_y(n));
    
    % Compute the normal for the brightest point of the sphere
    norm(1,n) = [x(n) - xc] / r ;
    norm(2,n) = -[y(n) - yc] / r ;
    norm(3,n) =  - sqrt ( 1 - norm(1,n)^2 - norm(2,n)^2 );
    dot_p = norm(1,n)*camera_dir(1) + norm(2,n)*camera_dir(2) + norm(3,n)*camera_dir(3);
    L(:,n) = 2*(dot_p)*(norm(:,n)) - camera_dir;
  end
    
      
  return;

