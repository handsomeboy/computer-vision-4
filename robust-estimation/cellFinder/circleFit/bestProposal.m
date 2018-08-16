function [circle] = bestProposal(circles, sigmaGM, normals, p)
% [circle] = bestProposal(circles, sigmaGM, normals, p)
% Chose the best circle from the proposed circles.
% Input
%  circles K x 3 matrix each row containing [x0(1), x0(2), r] for 
%          the center and radius of a proposed circle.
%  sigmaGM - the robust scale parameter 
%  normals - P x 2 edgel normal data
%  p       - P x 2 edgel position data.
% Output
%  circle  1 x 3  best circle params [x0(1), x0(2), r] from the proposals
    
% Get size of inputs
numCircles = size(circles,1);
numEdgels = size(p,1);

% Initialize matrices
circle = zeros(1,3);

% Return if no circle given
if numCircles == 0
    return;
end

% Return given circle if only one circle given
if numCircles == 1
    circle = circles(1,:);
    return;
end

% Initialize min value and lowest index
min = 100000;
lowest_index = 1;

% Loop through and do the first iteration of GM robust estimation for each
% The one with the lowest distance between theta_new and theta_old is
% chosen as the best estimator

for i = 1:numCircles
    % Find the parameters required to solve a/b old
    % Old parameters
    xc = circles(i,1:2)';
    r = circles(i,3);
    xc_sq = sum(xc.^2);
    xc_sq = repmat(xc_sq, 1, numEdgels);
    r = repmat(circles(i,3),1,numEdgels);
    
    % Solve for the a/b old parameters
    a = -2*xc;
    b = xc_sq - r.^2;
    b = b';
    
    % Solve for xk_sq
    xk_sq = sum((p .^2)')';
    
    % Solve for ek 
    ek = p * a + b + xk_sq;
    ek_sq = ek.^2;
    
    % Solve for W (weight) matrix    
    w = (sigmaGM./(sigmaGM^2 + ek_sq)).^2;
    W = diag(w);
    maxW = max(w);
    
    % Construct the X matrix    
    % Construct the X matrix
    X = ones(numEdgels,3);
    X(:, [1,2]) = p;
    
    % Solve for the model parameters: theta
    theta = (X' * W * X) \ (-X' * W * xk_sq );
    
    % Solve for a/b
    a = theta([1 2]);
    b = theta(3);
    
    % Solve for new xc, and r
    xc = a ./ (-2);
    xc_sq = sum(xc .^2);
    r = sqrt(xc_sq - b);
    
    % Check for convergence    
    theta_old_sum = sum(circles(i,:));
    theta_sum = sum(theta);
    
    % If the new sum is less than the old sum, keep the new sum   
    if abs ( theta_sum - theta_old_sum ) < min
        lowest_index = i;
        min = theta_sum;
    end
end

%Output
circle = circles(lowest_index,:);        
return;