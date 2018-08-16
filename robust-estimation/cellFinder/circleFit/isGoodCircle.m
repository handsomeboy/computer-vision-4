function [goodCircle] = isGoodCircle(x0, r, w, ...
                                     circleEstimates, nFound)
  % [goodCircle] = isGoodCircle(x0, r, w, normals, ...
  %                                  circleEstimates, nFound)
  % Decide if the circle with parameters x0 and r, with fitted
  % weights w, is to be added to the current set of circles.
  % Input:
  %  x0 2-vector for circle center
  %  r  radius of circle
  %  w  robust weights of edgels supporting this circle,
  %     weights have been scaled to have max possible value 1.
  %  circleEstimates C x 3 array of circle parameters, the first
  %     nFound rows of which specifies the [x0(1), x0(2) r] parameters
  %     for a previously fitted circle for this data set.
  %  nFound the number of previously fitted circles stored in circleEstimates
  % Output:
  %  goodCircle boolean, true if you wish to add this circle to the
  %             list of previously fitted circles.
  
  x0 = x0(:)';  % Decide, row or column
  
  w_Thr = 0.3*2*pi*r;
  d_Thr = 0.25*r;
  r_max = 0.6*r;
  
  % If no previous circles are found
  if nFound == 0
      % If the sum of the weights is greater than the threshold
      if sum(w) > w_Thr
          goodCircle = true;
      else
          goodCircle = false;
      end
  % Loop for if previous circles were found    
  else
      % Loop through number of circles found
      for i = 1:nFound
          % Calculate distance between centers
          d = pdist([x0;circleEstimates([1 2],i)']);
          
          % Find the circles that have a different in radii greater than
          % the distance between centers
          rMinusD(i,1) = (r + circleEstimates(3, i)) - d;
          
          % Store the radii of all the estimates
          rMinusD(i,2) = circleEstimates(3, i);
      end
      
      % Crop out only circles that interest 
      crop = rMinusD(:, 1) > 0;
      rMinusD = rMinusD(crop, :);
      
      % Crop out only circles that interest by less than 60%
      crop = rMinusD(:, 1) < r_max;
      rMinusD = rMinusD(crop, :);
      
      % Calculate the intersection radius
      rMinusD(:,3) = acos((r^2 + rMinusD(:,1).^2 - rMinusD(:,2).^2) ./ ...
                        (2 * rMinusD(:,1).* rMinusD(:,2)));
                    
      % Find the intersection angle between it and all the circles it
      % intersects with
      angle_diff = (2*pi) - 2*sum(rMinusD(:,3));
      
      % compare with threshold
      if (sum(w) > angle_diff*d_Thr)
          goodCircle = true;
      else
          goodCircle = false;
      end
  end
  
return;  


      
      
      
  