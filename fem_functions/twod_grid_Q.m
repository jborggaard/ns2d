function [ xy, bnodes ] = twod_grid_Q( jmax, kmax, jdim, kdim )
%TWOD_GRID_Q Generates a structured uniform quad grid.
%-------------------------------------------------------------------------------
%  Description - Creates a uniform grid on a quad domain.
%                
%
%  Version: 1.1
%
%  Usage:   [ x ] = twod_grid_Q( jmax, kmax, jdim, kdim )
%
%  Input:   (jmax,kmax)
%               Maximum value of the domain in each direction.  
%                 D = [0,jmax] x [0,kmax]
%           (jdim,kdim)
%               Number of nodes in each direction.
%
%  Output:
%           xy     - nx2 matrix with [x; y] coordinates. 
%           bnodes - n vector of boundary node markers with
%                      0 - internal node
%                      1 - edge node
%                      2 - corner node
%
%  Copyright (c) 2015, Alan Lattimer, Virginia Tech
%
%  Changes
%     June 2015 (1.1) - Optimized by removing for loops
%                     - Now outputs boundary node marker instead of
%                       internal node marker 
%-------------------------------------------------------------------------------
  
  tmp= repmat(linspace(0,jmax,jdim)',kdim,1);
  xy(:,1) = tmp(:);

  tmp= repmat(linspace(0,kmax,kdim),jdim,1);
  xy(:,2) = tmp(:);

  bnodes = sum([xy(:,1)==0, xy(:,1)==jmax, xy(:,2)==0, xy(:,2)==kmax],2);
  
end

