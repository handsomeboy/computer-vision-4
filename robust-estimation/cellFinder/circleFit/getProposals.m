function [circles] = getProposals(normals, p, numGuesses)
  % [circles] = getProposals(normals, p, numGuesses)
  % Attempt to produce up to numGuesses circle proposals from
  % the edgel data p and normals.  For typical data sets
  % we will be able to produce numGuesses proposals.  However,
  % on some datasets, say with only a few edgels, we may not be
  % able to generate any proposals.  In this case size(circles,1)
  % can be zero.
  % Input:
  %  normals - N x 2 edgel normals
  %  p         N x 2 edgel positions
  %  numGuesses - attempt to propose this number of circles.
  % Return:
  %   circles a P x 3 array, each row contains [x0(1) x0(2) r]
  %           with 0 <= P <= numGuesses.
 
 %Initialize parameters
 circles = zeros(0,3);
 
 % Get size of edgels given to function
 N = size(p,1);
 
 % Check if two few edgels are received
 if N < 2
     disp('Error: Not enough edgels');
     return;
 end;
 
 % Find x min, max, and range
 x_min = min(p(:,1));
 x_max = max(p(:,1));
 x_range = x_max - x_min;
 
 % Find y min, max, and range
 y_min = min(p(:,2));
 y_max = max(p(:,2));
 y_range = y_max - y_min;
 
 % Find size of edgel cluster
 edgels_size = x_range*y_range;
 
 % Control variables
 R_max = 0.5*max(x_range,y_range);
 R_min = 0.01*min(x_range,y_range);
 tolerance_ratio = 0.001;
 C_tol = tolerance_ratio*edgels_size;
 R_tol = tolerance_ratio*edgels_size;
 dist_min = 5;
 dist_max = 40;
 
 % Counter to enter found circles in order in circle matrix
 count = 1; 
 
 % Try to get proposals numGuesses times
 for n = 1:numGuesses

     % Randomly sample a point in the samples
     edgel_indx = randsample(1:N, 1);
     
     % Acquire location and normal of the first random index
     P1x = p(edgel_indx,1);
     P1y = p(edgel_indx,2);
     N1x = normals(edgel_indx,1);
     N1y = normals(edgel_indx,2);
       
     % Initialize variables for while loop
     distance = 0;
     next_indx = edgel_indx;
     P2y = 0;
     P2x = 0;
     
     % Counter to break loop if no good point found
     distance_check_count = 0;
     
     % Find point within distance tolerance
     
     while ( distance < dist_min || distance  > dist_max )
         next_indx = next_indx + 10;
         
         % Correct if index surpasses size of p
         if next_indx > N
             next_indx = next_indx - N;
         end
         
         % Acquire location of new index
         P2x = p(next_indx,1);
         P2y = p(next_indx,2);
         
         % Calculate distance between points
         distance = sqrt((P2x-P1x)^2+(P2y-P1y)^2);
         
         % Break if no good point found within 10 loops
         distance_check_count = distance_check_count + 1;
         if distance_check_count >= 10
             break;
         end
     end
     
     % Find normals for new points
     N2x = normals (next_indx,1);
     N2y = normals (next_indx,2);
       
     % Use x parameters to determine the radius
     rx = abs((P2x-P1x)/(N1x-N2x));
     
     % Use y parameters to determine the radius
     ry = abs((P2y-P1y)/(N1y-N2y));
     
     % Continue if found radii is greater than half the bounding window
     if rx > R_max || ry > R_max || rx < R_min || ry < R_min ...
            || abs(ry-rx) > R_tol
         continue;
     end
          
     % Find the average radius
     ravg = (rx+ry)/2;
     
     % Calculate the center of the circle using the first random index
     C1x = P1x + ravg*N1x;
     C1y = P1y + ravg*N1y;
     center1 = [C1x C1y];
     
     % Calculate the center of the circle using the second random index
     C2x = P2x + ravg*N2x;
     C2y = P2y + ravg*N2y;
     center2 = [C2x C2y];
     
     %Calculate the euclidean distance between the found centers
     C_distance = sqrt( (C2y-C1y)^2 + (C2x-C1x)^2 );
     
     % Break if distance between found centers are too far apart
     if C_distance > C_tol
         continue;
     end
     
     % Find the average circle
     Cx = (C1x+C2x)/2;
     Cy = (C1y+C2y)/2;
     center = [Cx Cy];
     
     % Add to circles
     circles(count,:) = [Cx, Cy, ravg];
     count = count + 1;
 end
 
end
     
  
  
