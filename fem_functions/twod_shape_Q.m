function [x_g,w_g,phi,p_x,p_y] = twod_shape_Q(x,r,s,w)
%-------------------------------------------------------------------------------
%  twod_shape_Q.m - computes test functions and derivatives on Q1,Q2,...
%                   elements given element coordinates.  Values are returned
%                   at the quadrature points (r,s)  assumed given in 
%                   (-1,1)x(-1,1)
%
%  Copyright (c) 2011, Jeff Borggaard, Virginia Tech
%  Version: 1.1
%
%  Usage:    [x_g,w_g,phi,p_x,p_y] = twod_shape_Q(x_local,r,w)
%
%  Variables:     x_local
%                         Coordinates of the element nodes
%                 (r,s)
%                         Coordinates of quadrature points in unit triangle
%                 w
%                         Quadrature weights associated with (r,s)
%
%                 x_g
%                         Coordinates of quadrature points in the element
%                 w_g
%                         Quadrature weights scaled by the element Jacobian
%                 phi
%                         Value of element shape functions at x_g
%                 p_x
%                 p_y
%                         First spatial derivatives of phi
%
% Modifications:
%  May 2015 - implemented quadratic elements
%-------------------------------------------------------------------------------
  %% set up
  [n,ncoord] = size(x);

%  fprintf('Space Dimension           : %i\n',ncoord)
%  fprintf('Number of shape functions : %i\n',n)

  if (ncoord ~= 2)
     fprintf('### WARNING : TWOD_SHAPE_Q ###\n')
     fprintf('Number of space coordinates : %i\n',ncoord)
     fprintf('### WARNING : TWOD_SHAPE_Q ###\n')
  end

  rule = length(r);
  phi = zeros(rule,n);  p_r = zeros(rule,n);  p_s = zeros(rule,n);
   
  %  Compute shape function and derivatives at quadrature points
  %% ---------------------------------------------------------------------------
  if (n == 4)
    %  4-noded bi-linear quadralateral element
    %%--------------------------------------------------------------------------
    phi(:,1) = (r-1).*(s-1)/4;  p_r(:,1) = (s-1)/4;  p_s(:,1) = (r-1)/4;
    phi(:,2) =-(r+1).*(s-1)/4;  p_r(:,2) =-(s-1)/4;  p_s(:,2) =-(r+1)/4;
    phi(:,3) = (r+1).*(s+1)/4;  p_r(:,3) = (s+1)/4;  p_s(:,3) = (r+1)/4;
    phi(:,4) =-(r-1).*(s+1)/4;  p_r(:,4) =-(s+1)/4;  p_s(:,4) =-(r-1)/4;
    
    x_g = phi*x;
    
    x_r = p_r*x;
    x_s = p_s*x;
    
    jac = x_r(:,1).*x_s(:,2) - x_s(:,1).*x_r(:,2);
    w_g = jac.*w;
   
    rx = x_s(:,2)./jac;
    sx =-x_r(:,2)./jac;
    ry =-x_s(:,1)./jac;
    sy = x_r(:,1)./jac;
    
    p_x = diag(rx)*p_r + diag(sx)*p_s;
    p_y = diag(ry)*p_r + diag(sy)*p_s;
    
  
  elseif (n == 8)
    %  8-noded serendipity element  (counterclockwise numbering)
    %
    %     4---7---3
    %     |       |   ( 1, x, x^2, xy, x^2y, xy^2, y, y^2 )
    %     8       6
    %     |       |
    %     1---5---2
    %
    %  hold off on 1/4 factor until the end.
    phi(:,1) =     r.*(1-r).*    s.*(1-s);   p_r(:,1) = (1-2*r).*s.*(1-s);
                                             p_s(:,1) = r.*(1-r).*(1-2*s);
    phi(:,2) = (1+r).*    r.*    s.*(1-s);   p_r(:,2) = (1+2*r).*s.*(1-s);
                                             p_s(:,2) = (1+r).*r.*(1-2*s);
    phi(:,3) = (1+r).*    r.*(1+s).*    s;   p_r(:,3) = (1+2*r).*(1+s).*s;
                                             p_s(:,3) = (1+r).*r.*(1+2*s);
    phi(:,4) =     r.*(1-r).*(1+s).*    s;   p_r(:,4) = (1-2*r).*(1+s).*s;
                                             p_s(:,4) = r.*(1-r).*(1+2*s);
    phi(:,5) = (1+r).*(1-r).*    s.*(1-s);   p_r(:,5) =    -2*r.*s.*(1-s);
                                             p_s(:,5) = (1-r.^2).*(1-2*s);
    phi(:,6) = (1+r).*    r.*(1+s).*(1-s);   p_r(:,6) = (1+2*r).*(1-s.^2);
                                             p_s(:,6) = (1+r).*r.* (-2*s);
    phi(:,7) = (1+r).*(1-r).*(1+s).*    s;   p_r(:,7) =    -2*r.*(1+s).*s;
                                             p_s(:,7) = (1-r.^2).*(1+2*s);
    phi(:,8) =     r.*(1-r).*(1+s).*(1-s);   p_r(:,8) = (1-2*r).*(1-s.^2);
                                             p_s(:,8) = r.*(1-r).* (-2*s);

    phi = phi/4;  p_r = p_r/4;  p_s = p_s/4;
    
    x_g = phi*x;
    
    x_r = p_r*x;
    x_s = p_s*x;
    
    jac = x_r(:,1).*x_s(:,2) - x_s(:,1).*x_r(:,2);
    w_g = jac.*w;
   
    rx = x_s(:,2)./jac;
    sx =-x_r(:,2)./jac;
    ry =-x_s(:,1)./jac;
    sy = x_r(:,1)./jac;
    
    p_x = diag(rx)*p_r + diag(sx)*p_s;
    p_y = diag(ry)*p_r + diag(sy)*p_s;
   
    
  elseif (n == 9)
%    error('element not implemented')
    %  9-noded quadratic (with (1-x^2)*(1-y^2) term)
    
    phi(:,1) = (r.^2 - r) .* (s.^2 - s) ./ 4;  
    phi(:,2) = (r.^2 + r) .* (s.^2 - s) ./ 4;
    phi(:,3) = (r.^2 + r) .* (s.^2 + s) ./ 4;
    phi(:,4) = (r.^2 - r) .* (s.^2 + s) ./ 4;
    phi(:,5) = (1 - r.^2) .* (s.^2 - s) ./ 2;
    phi(:,6) = (r.^2 + r) .* (1 - s.^2) ./ 2;
    phi(:,7) = (1 - r.^2) .* (s.^2 + s) ./ 2;
    phi(:,8) = (r.^2 - r) .* (1 - s.^2) ./ 2;
    phi(:,9) = (1 - r.^2) .* (1 - s.^2);
    
                                             
    p_r(:,1) = (2*r - 1) .* (s.^2 - s) ./ 4;
    p_r(:,2) = (2*r + 1) .* (s.^2 - s) ./ 4;
    p_r(:,3) = (2*r + 1) .* (s.^2 + s) ./ 4;
    p_r(:,4) = (2*r - 1) .* (s.^2 + s) ./ 4;
    p_r(:,5) = (   -2*r) .* (s.^2 - s) ./ 2;
    p_r(:,6) = (2*r + 1) .* (1 - s.^2) ./ 2;
    p_r(:,7) = (   -2*r) .* (s.^2 + s) ./ 2;
    p_r(:,8) = (2*r - 1) .* (1 - s.^2) ./ 2;
    p_r(:,9) = (   -2*r) .* (1 - s.^2);
    
    p_s(:,1) = (r.^2 - r) .* (2*s - 1) ./ 4;
    p_s(:,2) = (r.^2 + r) .* (2*s - 1) ./ 4;
    p_s(:,3) = (r.^2 + r) .* (2*s + 1) ./ 4;
    p_s(:,4) = (r.^2 - r) .* (2*s + 1) ./ 4;
    p_s(:,5) = (1 - r.^2) .* (2*s - 1) ./ 2;
    p_s(:,6) = (r.^2 + r) .* (   -2*s) ./ 2;
    p_s(:,7) = (1 - r.^2) .* (2*s + 1) ./ 2;
    p_s(:,8) = (r.^2 - r) .* (   -2*s) ./ 2;
    p_s(:,9) = (1 - r.^2) .* (   -2*s);


    x_g = phi*x;
    
    x_r = p_r*x;
    x_s = p_s*x;
    
    jac = x_r(:,1).*x_s(:,2) - x_s(:,1).*x_r(:,2);
    w_g = jac.*w;
   
    rx = x_s(:,2)./jac;
    sx =-x_r(:,2)./jac;
    ry =-x_s(:,1)./jac;
    sy = x_r(:,1)./jac;
    
    p_x = diag(rx)*p_r + diag(sx)*p_s;
    p_y = diag(ry)*p_r + diag(sy)*p_s;

  else
    fprintf('element not supported');
  end
