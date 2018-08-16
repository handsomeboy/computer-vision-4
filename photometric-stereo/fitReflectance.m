function [n, albedo] = fitReflectance(im, L)
  % [n, albedo] = fitReflectance(im, L)

    [nPix,nDirChrome] = size(im);
    
    g = (im*L')/(L*L');
    n = normr(g);
    albedo = sqrt(g(:,1).^2 + g(:,2).^2 + g(:,3).^2);

  % Input:
  %   im - nPix x nDirChrome array of brightnesses,
  %   L  - 3 x nDirChrome array of light source directions.
  % Output:
  %   n - nPix x 3 array of surface normals, with n(k,1:3) = (nx, ny, nz)
  %       at the k-th pixel.
  %   albedo - nPix x 1 array of estimated albdedos
    
  return;


