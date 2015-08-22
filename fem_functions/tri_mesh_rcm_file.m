function [] = tri_mesh_rcm(mesh_root)
%%
%  tri_mesh_rcm - a function that reads in a triangular_mesh and produces a
%                 renumbering strategy based on the reverse Cuthill-McKee
%                 algorithm.  The tri_mesh can correspond to either 
%                 linear (3 vertices with right-hand orientation) or 
%                 quadratic (3 vertices + 3 midside nodes).
%
%%

  x       = load(strcat(mesh_root,'.node'));
  e_conn  = load(strcat(mesh_root,'.elem'));

  [n_node,n_dim] = size(x);
  [n_elem,n_dof] = size(e_conn);


  %  preallocate storage for the adjacency graph
  
  fprintf('building adjacency graph\n')
%  adjacency = sparse(n_node,n_node);

  II = zeros(n_dof*n_dof*n_elem,1);
  JJ = II;
  AA = ones(n_dof*n_dof*n_elem,1);
  ntriple = 0;
  
  for n_el = 1:n_elem
    for i=1:n_dof
      for j=1:n_dof
        ntriple = ntriple + 1;
        II(ntriple) = e_conn(n_el,i);
        JJ(ntriple) = e_conn(n_el,j);
      end
    end
  end

  adjacency = sparse(II,JJ,AA,n_node,n_node);
  
  figure(12)
  spy(adjacency)

  fprintf('calling symrcm\n')

  p = symrcm(adjacency);

  for n=1:n_node
    pinv(p(n)) = n;
  end

  figure(13)
  spy(adjacency(p,p))

  clear adjacency

  x_new = x(p,:);
  clear x

  for n_el = 1:n_elem
    e_conn_new(n_el,:) = pinv(e_conn(n_el,:));
  end


  mesh_root_new = strcat(mesh_root,'rcm');
 
  write_mesh_2d(x_new,e_conn_new,mesh_root_new)


end % function  tet_mesh_boundary
  
