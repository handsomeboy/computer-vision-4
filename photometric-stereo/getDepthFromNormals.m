function [depth] = getDepthFromNormals(n, mask)
%    depth = getDepthFromNormals(n,mask);

    [N,M] = size(mask);
     
    % Initialize arrays
    idx = zeros(N,M);
    tmp = zeros(N,M);
    depth = zeros(N,M);
    
    % Find values within the mask that have values
    [ row_num, col_num ] = find ( mask );
    numPixels = size ( row_num, 1 );
    [ elements_w_values ] = zeros ( numPixels, 2 ); % Combine above into one matrix
    elements_w_values = [ row_num col_num ];
    
    % Create sparse matrices
    A = sparse(2*numPixels, numPixels);
    v = sparse(2*numPixels, 1);
    
    % Create an index of elements to manipulate with n values greater than
    % 0
    
    for i = 1:numPixels
        idx ( elements_w_values(i,1), elements_w_values(i,2)) = i;
    end
    
    row_wo_neighbours = [];
    for i = 1:numPixels
        
        % Get pixel element location
        I = elements_w_values (i, 1);
        J = elements_w_values (i, 2);
        
        % Normal values at that pixel element
        Nx = n(I, J, 1);
        Ny = n(I, J, 2);
        Nz = n(I, J, 3);
        
        % Index locations for A and v matrices
        v_idx = 2*i - 1; % Odd
        h_idx = 2*i; % Even
        
   % First-order Forward Differences
     % Vertical Calculations
        % Check to see to element below exists
        if ( idx(I+1,J) > 0 )
            adj_index = idx(I+1, J); % Adjacent neighbour's index number
            
            A(v_idx, i) = -Nz;
            A(v_idx, adj_index ) = Nz;
            v(v_idx) = Ny;
        end   
                    
     % Horizontal Calculations
        % Check if element to the right exists
        if ( idx(I,J+1) > 0 )
            adj_index = idx(I,J+1); % Adjacent neighbour's index number
            
            A(h_idx, i) = - Nz;
            A(h_idx, adj_index) = Nz;
            v(h_idx) = - Nx;
        end
    end
    
    % Calculate z
    z = A \ v;
    
    % Rescale z from 0 to 255
    z = ( z - min(z) ) / ( max(z) - min(z) ) * 255;
    
    % Invert z so closest point has a depth of zero and furthest point
    % has a depth of 255
    z = 255 - z;
       
    for i = 1:numPixels
        I = elements_w_values (i, 1);
        J = elements_w_values (i, 2);
        depth(I,J) = z(i,1);
    end
    
  % Input:
  %    n is an [N, M, 3] matrix of surface normals (or zeros
  %      for no available normal).
  %    mask logical [N,M] matrix which is true for pixels
  %      at which the object is present.
  % Output:
  %    depth an [N,M] matrix providing depths which are
  %          orthogonal to the normals n (in the least
  %          squares sense).
  %

end