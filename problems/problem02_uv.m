function [x,e_conn,ide_u,ide_p,dir_u] = problem2_uv(N1,N2)

  etype = 'quadratic';
  [x,e_conn,b_node] = twod_mesh(0,1,0,5,etype,N1,N2);

  [n_nodes   , n_dimensions] = size(x     );
  [n_elements, nel_dof     ] = size(e_conn);
  

  inflow_counter  = 0;
  outflow_counter = 0;
  wall_counter    = 0;
    
  for n = b_node
      
    % Inflow Boundary Nodes
    if ( abs(x(n,2)-5) < 1e-8 )
      inflow_counter = inflow_counter + 1;
      inflow_nodes(inflow_counter) = n;
      
    % Outflow Boundary Nodes
    elseif ( abs(x(n,2)-0) < 1e-8 && ( x(n,1)>0 && x(n,1)<1 ) )
      outflow_counter = outflow_counter + 1;
      outflow_nodes(outflow_counter) = n;
      
    % Wall Boundary Nodes
    else
      wall_counter = wall_counter + 1;
      wall_nodes(wall_counter) = n;       
    end
      
  end

  n_diru      = 0;
  n_equations = 0;
  ide_u       = zeros(n_nodes,2);
  ide_p       = zeros(n_nodes,1);

  dir_u       = zeros(2*length(b_node),1);
  
  for n_nd=inflow_nodes
    n_diru        = n_diru + 1;
    ide_u(n_nd,1) = -n_diru;
    dir_u(n_diru) = 0;

    n_diru        = n_diru + 1;
    ide_u(n_nd,2) = -n_diru;
    dir_u(n_diru) = - ( x(n_nd,1) )*( 1-x(n_nd,1) ) / 0.25;
  end
 
  %  One option for the outflow boundary condition, another
  %  option would be to comment out the "for loop" below.
  for n_nd=outflow_nodes
    n_diru        = n_diru + 1;
    ide_u(n_nd,1) = -n_diru;
    dir_u(n_diru) = 0;
  end
  
  
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
  
  dir_u = dir_u(1:n_diru);
  
end
