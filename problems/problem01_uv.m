function [x,e_conn,ide_u,ide_p,dir_u] = problem01_uv(N1,N2)
%----------------------------------------------------------------------------78-
%  problem1_uv - a function that provides input data to a 2D Navier-Stokes
%                finite element solver (unknown variables u, v, and p).  
%
%                this simulates the Poiseuille flow in a 5x1 channel
%                (flow from left to right)
%
%  useage:
%        [x,e_conn,ide_u,ide_p,dir_u] = problem01_uv(N1,N2)
%
%  inputs: 
%        N1, N2
%                numbers of points in the horizontal and vertical directions
%                 
%-------------------------------------------------------------------------------


  etype = 'quadratic';
  [x,e_conn,b_node] = twod_mesh(0,5,0,1,etype,N1,N2);
    
  
  inflow_counter  = 0;
  outflow_counter = 0;
  wall_counter    = 0;
    
  for n = b_node
      
    % Inflow Boundary Conditions
    if ( abs(x(n,1)-0) < 1e-8 )
      inflow_counter = inflow_counter + 1;
      inflow_nodes(inflow_counter) = n;
      
    % Outflow Boundary Conditions
    elseif ( abs(x(n,1)-5) < 1e-8 && ( x(n,2)>0 && x(n,2)<1 ) )
      outflow_counter = outflow_counter + 1;
      outflow_nodes(outflow_counter) = n;
      
    % Wall Boundary Conditions
    else
      wall_counter = wall_counter + 1;
      wall_nodes(wall_counter) = n;       
    end
      
  end


  [n_nodes   , n_dimensions] = size(x     );
  [n_elements, nel_dof     ] = size(e_conn);
  
  %-----------------------------------------------------------------------------
  %  Determine equation numbers, set up boundary condition information
  %-----------------------------------------------------------------------------
  n_diru      = 0;
  n_equations = 0;
  ide_u       = zeros(n_nodes,2);
  ide_p       = zeros(n_nodes,1);

  dir_u       = zeros(2*length(b_node),1);
    
  %  Quadratic (fully-developed) inflow profile 
  for n_nd=inflow_nodes
    n_diru        = n_diru + 1;
    ide_u(n_nd,1) = -n_diru;
    dir_u(n_diru) = 4*( x(n_nd,2)-x(n_nd,2)^2 );

    n_diru        = n_diru + 1;
    ide_u(n_nd,2) = -n_diru;
    dir_u(n_diru) = 0;
  end

  %  One option for the outflow boundary condition, another
  %  option would be to comment out the "for loop" below.
  for n_nd=outflow_nodes
    n_diru        = n_diru + 1;
    ide_u(n_nd,2) = -n_diru;
    dir_u(n_diru) = 0;
  end
  
  %  No slip, no penetration boundary conditions on the walls
  for n_nd=wall_nodes
    n_diru        = n_diru + 1;
    ide_u(n_nd,1) = -n_diru;
    dir_u(n_diru) = 0;

    n_diru        = n_diru + 1;
    ide_u(n_nd,2) = -n_diru;
    dir_u(n_diru) = 0;
  end

  % loop over all nodes, determine if dirichlet boundary condition number 
  % has been set.  If not, it is a degree of freedom, so number it.
  for n_nd=1:n_nodes
    if ( ide_u(n_nd,1)==0 )
      n_equations   = n_equations + 1;
      ide_u(n_nd,1) = n_equations;
    end
    
    if ( ide_u(n_nd,2)==0 )
      n_equations   = n_equations + 1;
      ide_u(n_nd,2) = n_equations;
    end
  end

  for n_el=1:n_elements
    vertex = e_conn(n_el,1:3);
    if (ide_p(vertex(1))==0)
      n_equations      = n_equations + 1;
      ide_p(vertex(1)) = n_equations;
    end
    if (ide_p(vertex(2))==0)
      n_equations      = n_equations + 1;
      ide_p(vertex(2)) = n_equations;
    end
    if (ide_p(vertex(3))==0)
      n_equations      = n_equations + 1;
      ide_p(vertex(3)) = n_equations;
    end
  end
  
  %  Compress the dirichlet boundary condition array
  dir_u = dir_u(1:n_diru);
  
end
