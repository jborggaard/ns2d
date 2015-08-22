function [ x, e_conn, var, v_indices ] = twod_quad_to_lin( x, e_conn, var )
%TWOD_QUAD_TO_LIN Reduces a quadratic mesh to a linear mesh.
%   Input geometry (and optional variables) for a quadratic tetrahedral
%   mesh.  This routine finds the vertices and returns the geometry
%   (and optional variables) for the corresponding linear tetrahedral mesh.

%  Make a linear mesh from a quadratic mesh 
  [n_nodes   , n_dimensions] = size(x     );
  [n_elements, nel_dof     ] = size(e_conn);
  
  is_vertex = zeros(n_nodes);
  for n=1:n_elements
    for d=1:3
      is_vertex(e_conn(n,d)) = 1;
    end
  end
  
  % gets a list of the quadratic nodes to be kept in linear mesh
  v_indices = find(is_vertex);  
 
  l_count = 0;
  for n=1:n_nodes
    if ( is_vertex(n) )
      l_count = l_count + 1;
      is_vertex(n) = l_count;
    end
  end
  
  x = x(v_indices,:);
  e_conn = is_vertex(e_conn(:,1:3));

  var = var(v_indices,:);
end

