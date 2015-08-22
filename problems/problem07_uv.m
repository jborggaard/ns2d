function [x,e_conn,ide_u,ide_p,dir_u] = problem07_uv(N,ifsens)
%% -----------------------------------------------------------------------------
%  Time dependent verification example for a 2d Navier-Stokes solver.  This
%  solution represents the Stommel ocean circulation model.  A time dependent
%  solution and corresponding manufactured forcing function can be found in
%  the files  exact_solution7.m  and  f_function7_2d.m, respectively.
%%------------------------------------------------------------------------------
  if ( nargin<2 )
    ifsens = 0;
  end

  R      =   1.5198e-6;     % nondimensional friction coefficient
  beta   =   0         ;    % nondimensional Coriolis effect parameter
  x_len  =   1.591549430918954;
  y_len  =   1;

  [x,e_conn]   = twod_mesh(0,x_len,0,y_len,'quadratic',1.5*(N-1)+1,N);
  
  [x,e_conn]   = tri_mesh_rcm(x,e_conn);  
  [e_conn,adj] = tri_mesh_corner(e_conn);
  [b_node]     = tri_mesh_boundary(e_conn,adj);
  
  exact_counter    = 0;
  exact_nodes      = zeros(size(b_node));
  for n = b_node
            
    % Exact Boundary Conditions
    if ( abs(x(n,1)-0) < 1e-8 || abs(x(n,1)-x_len) < 1e-8 || ...
         abs(x(n,2)-0) < 1e-8 || abs(x(n,2)-y_len) < 1e-8 )
      exact_counter = exact_counter + 1;
      exact_nodes(exact_counter) = n;       
    end
      
  end
  exact_nodes = exact_nodes(1:exact_counter);
    
  n_nodes    = size(x     ,1);
  n_elements = size(e_conn,1);
  
  %-----------------------------------------------------------------------------
  %  Determine equation numbers, set up boundary condition information
  %-----------------------------------------------------------------------------
  n_diru      = 0;
  n_equations = 0;
  ide_u       = zeros(n_nodes,2);
  ide_p       = zeros(n_nodes,1);
  time        = 0;

  if ( ifsens )
    dR     = 1e-6*R;
    dbeta  = 1e-6;
 
    dir_u     = zeros(2*length(exact_nodes),3);
    
    for n_nd=exact_nodes
      [u  ,v] = exact_solution7(x(n_nd,:),time,R,beta);
      [upR, vpR] = exact_solution7(x(n_nd,:),time,R+dR,beta);
      [upb, vpb] = exact_solution7(x(n_nd,:),time,R,beta+dbeta);
      
      n_diru          = n_diru + 1;
      ide_u(n_nd,1)   = -n_diru;
      dir_u(n_diru,1) = u;
      dir_u(n_diru,2) = (upR-u)/dR;
      dir_u(n_diru,3) = (upb-u)/dbeta;

      n_diru          = n_diru + 1;
      ide_u(n_nd,2)   = -n_diru;
      dir_u(n_diru,1) = v;
      dir_u(n_diru,2) = (vpR-v)/dR;
      dir_u(n_diru,3) = (vpb-v)/dbeta;
    end
    
  else
    dir_u     = zeros(2*length(exact_nodes),1);
  
    for n_nd=exact_nodes
      
      [u, v]          = exact_solution7(x(n_nd,:),time);
      n_diru          = n_diru + 1;
      ide_u(n_nd,1)   = -n_diru;
      dir_u(n_diru,1) = u;

      n_diru          = n_diru + 1;
      ide_u(n_nd,2)   = -n_diru;
      dir_u(n_diru,1) = v;
    end
  end
  
  
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
  
end % function
