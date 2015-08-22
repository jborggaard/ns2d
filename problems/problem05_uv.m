function [x,e_conn,ide_u,ide_p,dir_u] = problem05_uv(mesh_root)

% Flow over a Reuleaux shape (constant blockage under rotation).

    x      = load(strcat(mesh_root,'.node'));
    e_conn = load(strcat(mesh_root,'.elem'));
    
    b_node = load(strcat(mesh_root,'.boundary'));
    
    inflow_counter  = 0;
    outflow_counter = 0;
    wall_counter    = 0;
    
  % some pre-allocation
  inflow_nodes    = zeros(length(b_node),1);
  outflow_nodes   = inflow_nodes;
  wall_nodes      = inflow_nodes;
    
  for n = b_node'
        
    % Inflow Boundary Nodes
    if ( abs( x(n,1)+5 ) < 1e-8 )
      inflow_counter = inflow_counter + 1;
      inflow_nodes(inflow_counter) = n;

    % Outflow Boundary Nodes
    elseif ( abs( x(n,1)-15 ) < 1e-8 && x(n,2) > -5 + 1e-8 && x(n,2) < 5 - 1e-8 )
      outflow_counter = outflow_counter + 1;
      outflow_nodes(outflow_counter) = n;
        
    % Wall Boundary Nodes
    else
      wall_counter = wall_counter + 1;
      wall_nodes(wall_counter) = n;

    end
    
  end % boundary node list

  inflow_nodes = inflow_nodes(1:inflow_counter)';
  outflow_nodes = outflow_nodes(1:outflow_counter)';
  wall_nodes = wall_nodes(1:wall_counter)';
  
  
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
  

  for n_nd=inflow_nodes
    n_diru        = n_diru + 1;
    ide_u(n_nd,1) = -n_diru;
    dir_u(n_diru) = (5+x(n_nd,2)) * (5-x(n_nd,2)) / 25;

    n_diru        = n_diru + 1;
    ide_u(n_nd,2) = -n_diru;
    dir_u(n_diru) = 0 ;
  end

  for n_nd=wall_nodes
    n_diru        = n_diru + 1;
    ide_u(n_nd,1) = -n_diru;
    dir_u(n_diru) = 0;

    n_diru        = n_diru + 1;
    ide_u(n_nd,2) = -n_diru;
    dir_u(n_diru) = 0;
  end

  dir_u = dir_u(1:n_diru);

  
  % Number all of the unknowns
  
  numbering = 1;
  if (numbering)  % traditional numbering scheme (separate u and p)
    % loop over all nodes, determine if dirichlet boundary condition number 
    % has been set.  If not, it is a degree of freedom, so number it.
    for n_nd=1:n_nodes
      if ( ide_u(n_nd,1)==0 )
        n_equations   = n_equations + 1;
        ide_u(n_nd,1) = n_equations;

        n_equations   = n_equations + 1;
        ide_u(n_nd,2) = n_equations;
      end
    end

    for n_el=1:n_elements
      vertex = e_conn(n_el,1:3);
      for j=1:3
      if (ide_p(vertex(j))==0)
        n_equations      = n_equations + 1;
        ide_p(vertex(j)) = n_equations;
      end
      end
    end
    
  else
  % ALTERNATE NUMBERING (reduces entire bandwidth of Jacobian matrix)
    is_vertex = zeros(1,n_nodes);
    for n_el=1:n_elements
      vertex = e_conn(n_el,1:3);
      is_vertex(vertex) = [1 1 1];
    end
    
    for n_nd=1:n_nodes
      if ( ide_u(n_nd,1)==0 )
        n_equations   = n_equations + 1;
        ide_u(n_nd,1) = n_equations;
          
        n_equations   = n_equations + 1;
        ide_u(n_nd,2) = n_equations;
      end
      
      if ( is_vertex(n_nd)==1 )
%        if ( ide_p(n_nd,1)==0 )
          n_equations   = n_equations + 1;
          ide_p(n_nd,1) = n_equations;
%        end
      end
    end
    
  end
        
end % function
