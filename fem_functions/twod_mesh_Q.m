function [ x, e_conn, bnodes ] = twod_mesh_Q( order, jdim, kdim, domain )
%TWOD_MESH_Q Generates element connection matrix for a structured quad grid.
%-------------------------------------------------------------------------------
%  Description - Creates a connection matrix for a uniform mesh of quad
%                elements a regular rectangular domain.  Either linear 
%                quadratic quadrilateral elements are constructed
%
%  Version: 1.0
%
%  Usage:   [x, e_conn] = twod_mesh_Q(order, jmax,kmax,,domain)
%
%  Input:   order
%             order of element:
%                   order = 1 => linear quads
%                   order = 2 => quadratic quads
%
%           (jdim,kdim)
%             Number of nodes in each direction. 
%                 Requirements:
%                   Linear: at least 2 in each direction
%                   Quadratic: odd and at least 3 in each direction     
%                 Default: (11,11)
%
%           domain 
%              [jmax,kmax] - The domain is assumed to be
%                                  D = [0,jmax] x [0,kmax]
%              [jmin,jmax;min,ymax]
%
%  Output:
%           x 
%             nx2 matrix of the (j,k) coordinates of the quadrilateral grid
%
%           e_conn
%             Element connectivity matrix where each row is a
%             node on the quadrilateral element.
%
%           bnodes
%             List of the boundary nodes where
%               0 - internal node
%               1 - edge node
%               2 - corner node
%
%  Copyright (c) 2015, Alan Lattimer, Virginia Tech
%-------------------------------------------------------------------------------
  if (nargin<4)
    domain = [0,1;0,1];
  else
    [md,nd] = size(domain);
    if md == 1 && nd == 2
      domain = [[0;0] domain'];
    elseif md == 2 && md == 1
      domain = [[0;0] domain];
    elseif md ~= 2 || nd ~= 2
      error('Incorrect domain definition. Should be [jmax,kmax] or [jmin,jmax;kmin,kmax]')
    end
  end
  if (nargin < 3)
    warning('No dimensions set.  Defaulting to [11,11]');
    jdim = 11;
    kdim = 11;
  end

  xd = domain(1,:);
  yd = domain(2,:);

  [x, bnodes] = twod_grid_Q(1,1,jdim,kdim);

  x = [(x(:,1)*(xd(2)-xd(1)))+xd(1) (x(:,2)*(yd(2)-yd(1)))+yd(1)];

  if (order == 1)
    if (jdim<2 || kdim<2)
      error('Number of nodes must be >1 in each direction for linear elements.');
    end

    num_elem   = (jdim-1)*(kdim-1);
    % Each element has 4 nodes
    %   4------3
    %   |      |
    %   |      |
    %   1------2
    %
    e_conn = zeros(num_elem,4);

    dk2 = jdim;

    elem      = 1;
    for k = 1:kdim-1
      for j = 1:jdim-1       
        n1 = j + ((k-1)*dk2);
        n2 = n1 + 1;
        n3 = n2 + dk2;
        n4 = n1 + dk2;
        
        e_conn(elem,:) = [n1 n2 n3 n4];

        elem = elem + 1;
      end
    end
  elseif (order == 2)
    if (mod(jdim,2) == 0 || mod(kdim,2) == 0)
      error('Number of elements must be odd in every direction for quadratic elements.')
    end
    if (jdim<3 || kdim<3)
      error('Number of nodes must be at least 3 in each direction for quadratic elements.');
    end
    %    4------7------3  
    %    |             |
    %    |             |
    %    8      9      6
    %    |             |
    %    |             |
    %    1------5------2
    %
    %
    num_elem = (jdim-1)*(kdim-1)/4;
    e_conn = zeros(num_elem,9);

    dk1 = jdim;
    dk2 = 2*dk1;

    elem = 1;
    for k = 1:2:kdim-2
      for j = 1:2:jdim-2
        % Vertex Nodes
        n1 = j + ((k-1)*dk1);
        n2 = n1 + 2;
        n3 = n2 + dk2;
        n4 = n1 + dk2;
        e_conn(elem,1:4) = [n1 n2 n3 n4];

        % Edge Nodes
        n5 = n1 + 1;
        n6 = n2 + dk1;
        n7 = n4 + 1;
        n8 = n1 + dk1;

        e_conn(elem,5:8) = [n5 n6 n7 n8];

        % Center Node
        e_conn(elem,9) = n1+kdim+1;

        elem = elem + 1;
      end
    end
  else
    error('threed_hex_mesh:orderck','Order must be 1 or 2.  Only linear and quadratics implemented.');
  end






end % end threed_mesh_Q

