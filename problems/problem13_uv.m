function [x,e_conn,ide_u,ide_p,dir_u] = problem13_uv(override,msize)

  addpath('../fem_functions')
  
  if ( nargin==0 )
    override=false;
    msize   = 111814;
  end
  
  matlab_file = ['solutions_mat/ns_2d_13mesh_',sprintf('%06d',msize),'.mat'];

  verbose = true;    % flag for enabling command window output

  if ( exist(matlab_file,'file') && override==false )
    if ( verbose )
      fprintf(' loading solution from %s\n',matlab_file)
    end
    
    load (matlab_file, 'x', 'e_conn', 'boundary');

  else
    mesh_file = ['../../Meshes/cylinders_',sprintf('%06d',msize),'.msh'];
%    mesh_file = ['../../Meshes/cylinders_g1_25_',sprintf('%06d',msize),'.msh'];
      
    if ( verbose )
      fprintf(' importing mesh from %s',mesh_file)
    end
    
    [x,e_conn,boundary] = read_msh(mesh_file);

    if ( verbose )
      fprintf(' renumbering unknowns using symrcm\n')
    end
    
    [x,e_conn,p_inv] = tri_mesh_rcm(x,e_conn);
    for n=1:length(boundary)
      boundary(n).nodes = p_inv(boundary(n).nodes);
    end
    
    save(matlab_file,'x','e_conn','boundary');
  end
  
  n_nodes    = size(x     ,1);
  x = x(:,1:2);
  [n_elements,nel_dof] = size(e_conn);
  
  if ( nel_dof~=6 )
    error('problem13_uv: Expecting 2D Mesh with Quadratic Elements')
  end
  
%   figure; hold on
%   for n_nd=boundary(1).nodes
%     plot(x(n_nd,1),x(n_nd,2),'r+')
%   end
%   
%   for n_nd=boundary(2).nodes
%     plot(x(n_nd,1),x(n_nd,2),'g+')
%   end
%   
%   for n_nd=boundary(3).nodes
%     plot(x(n_nd,1),x(n_nd,2),'k*')
%   end
%   
%   for n_nd=boundary(4).nodes
%     plot(x(n_nd,1),x(n_nd,2),'b*')
%   end
  
  %-----------------------------------------------------------------------------
  %  Determine equation numbers, set up boundary condition information
  %-----------------------------------------------------------------------------
  %  Specify Dirichlet boundary conditions on velocities
  %    variable(1) = u
  %    variable(2) = v
  %    variable(3) = w
  
  %  Preallocation
  n_diru      = 0;
  n_equations = 0;
  ide_u       = zeros(n_nodes,2);
  ide_p       = zeros(n_nodes,1);

  dir_u       = zeros(length(boundary(1).nodes)...
                     +length(boundary(2).nodes)...
                     +length(boundary(3).nodes),1);

  %  Set inflow boundary conditions
  for n_nd=boundary(3).nodes
    n_diru        = n_diru + 1;
    ide_u(n_nd,1) = -n_diru;
    dir_u(n_diru) = 1;

    n_diru        = n_diru + 1;
    ide_u(n_nd,2) = -n_diru;
    dir_u(n_diru) = 0;
  end

  %  Set conditions on bottom cylinder
  for n_nd=boundary(1).nodes
    n_diru        = n_diru + 1;
    ide_u(n_nd,1) = -n_diru;
    dir_u(n_diru) = 0;

    n_diru        = n_diru + 1;
    ide_u(n_nd,2) = -n_diru;
    dir_u(n_diru) = 0;
  end

  %  Set conditions on top cylinder
  for n_nd=boundary(2).nodes
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

end % function
