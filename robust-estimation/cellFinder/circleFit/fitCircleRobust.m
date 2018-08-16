function [x0, rc, w, maxW] = fitCircleRobust(pts, initx0, initr, normals, sigmaGM)
%
% function [x0, r, w, maxW] = fitCircleRobust(pts, initx0, initr,
%                                  normals, sigmaGM)
%
%  minimize sum_i  rho( a'x_i + b + x_i'x_i, sigma)
%  w.r.t. a and b, where a = -2x_0 and b = x_0'x_0 - r^2
% Input:
%  pts: an Nx2 matrix containing 2D x_i, i=1...N
%  initx0: 2-vector initial guess for centre
%  initr: initial guess for radius 
%  normals: an Nx2 matrix of corresponding 2D edge normal directions nrml_i  
% Output
%  x0 - 2x1 fitted circle's center position
%  r  - fitted circle's radius
%  w  - N x 1 robust weights for edgels.
%  maxW  - maximum possible robust weight (i.e. for zero error)
%          This may be used by the calling code to scale the weights 
%          to be between 0 and 1.
  
% Initialize variables
x0 = initx0(:)';  % Make it a row vector.
xc = initx0;
rc = initr;
xk = pts;
theta_sum_old = 0;

% Find size of N
N = size(pts,1);

% Iteration Loop
for i=1:50
    
    % Find the parameters required to solve for a/b old
    xc_sq = sum(xc.^2);
    xc_sq = repmat(xc_sq, 1, N);
    rc = repmat(rc,1,N);
    
    % Solve for the old a/b parameter    
    a = -2*xc;
    b = xc_sq - rc.^2;
    b = b';
    
    % Solve for xk_sq
    xk_sq = sum((xk.^2)')';
    
    % Solve for ek     
    ek = xk*a + b + xk_sq;
    ek_sq = ek.^2;
    
    % Solve for the W(weight) matrix    
    w = (sigmaGM./(sigmaGM^2 + ek_sq)).^2;
    W = diag(w); % Construct the W matrix
    maxW = max(w); % find max weight value to return
    
    % Construct the X matrix
    X = ones(N,3);
    X(:, [1, 2]) = xk;
    
    % Solve for the new model parameters
    theta = (X'*W*X)\(-X'*W*xk_sq);
    
    % Solve for a/b from theta
    a = theta([1 2]);
    b = theta(3);
    
    % Solve for the new xc, and r
    xc = a ./ (-2);
    xc_sq = sum(xc.^2);
    rc = sqrt(xc_sq - b);
    
    %Check for convergence
    theta_sum = sum(theta);    
    if abs(theta_sum - theta_sum_old) < 0.05
        break;
    end
    theta_sum_old = theta_sum; % Set old values = new values
end

% Set values of the found circle
x0 = xc;
r = rc;

return;
    
    
    


