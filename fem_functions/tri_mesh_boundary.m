function [boundary_nodes,boundary_faces] = tri_mesh_boundary(e_conn,tri_adj)
%%------------------------------------------------------------------------------
%  tri_mesh_boundary - a function that reads in a tri_mesh along with neighbor
%                      information and produces lists of boundary nodes and 
%                      boundary edges.  The tri_mesh can correspond to either 
%                      linear (3 vertices with right-hand orientation) or 
%                      quadratic (3 vertices + 3 midside nodes).  In the
%                      latter case, we assume the following ordering:
%                        1  2  3  (1+2) (2+3) (1+3).
%
%  Usage:
%  [boundary_nodes,boundary_faces] = tri_mesh_boundary(e_conn,tri_adj)
%
%  Variables:
%            e_conn
%            tri_adj
%
%            boundary_nodes
%% -----------------------------------------------------------------------------

  [ boundary_faces ] = define_boundary_faces( e_conn, tri_adj );


  %  preallocate storage for boundary_nodes array
  [ n_faces, npf ] = size( boundary_faces );

  
  boundary_nodes = zeros( 1,n_faces*npf );
  boundary_nodes_counter = 0;

  for n_f=1:n_faces
    boundary_nodes( boundary_nodes_counter+1:boundary_nodes_counter+npf ) =...
      boundary_faces( n_f, 1:npf );
    boundary_nodes_counter = boundary_nodes_counter+npf;
  end

  [ boundary_nodes ] = unique( boundary_nodes );

end % function  tet_mesh_boundary
  



function [ boundary_faces ] = define_boundary_faces( e_conn, tri_adj )
%%
%  define_boundary_faces - a function that loops over elements, finds
%                          boundary faces from -1 entries to the
%                          adjacency array and stores information
%--------------------------------------------------------------------------
  [n_elem,n_dof] = size(e_conn);

  %  Estimate the number of boundary faces for preallocation
  %  (these will always be unique)
  if ( n_dof==3 )
    boundary_faces = zeros(3*n_elem,2);
  else 
    boundary_faces = zeros(3*n_elem,3);
  end
  
  boundary_faces_counter = 0;

  if ( n_dof==3 )
    %----------------------------------------------------------------------
    %  Case: Linear Triangle (3 vertices)
    %----------------------------------------------------------------------
    for n_el=1:n_elem
      % A 1-FACE
      if ( tri_adj(n_el,1) < 0 )
        boundary_faces_counter = boundary_faces_counter+1;
        boundary_faces(boundary_faces_counter,:) = ...
         [ e_conn(n_el, 2)  e_conn(n_el, 3) ];
      end


      % A 2-FACE
      if ( tri_adj(n_el,2) < 0 )
        boundary_faces_counter = boundary_faces_counter+1;
        boundary_faces(boundary_faces_counter,:) = ...
         [ e_conn(n_el, 3)  e_conn(n_el, 1) ];
      end


      % A 3-FACE
      if ( tri_adj(n_el,3) < 0 )
        boundary_faces_counter = boundary_faces_counter+1;
        boundary_faces(boundary_faces_counter,:) = ...
         [ e_conn(n_el, 1)  e_conn(n_el, 2) ];
      end

    end % element loop

    % end linear triangle case

  elseif ( n_dof==6 || n_dof==7 )
    %----------------------------------------------------------------------
    %  Case: Quadratic Triangle (6 noded or 7 noded Crouzeux-Raviart)
    %----------------------------------------------------------------------
    for n_el=1:n_elem
      % A 1-FACE
      if ( tri_adj(n_el,1) < 0 )
        boundary_faces_counter = boundary_faces_counter+1;
        boundary_faces(boundary_faces_counter,:) = ...
         [ e_conn(n_el, 2)  e_conn(n_el, 5)  e_conn(n_el, 3) ];
      end


      % A 2-FACE
      if ( tri_adj(n_el,2) < 0 )
        boundary_faces_counter = boundary_faces_counter+1;
        boundary_faces(boundary_faces_counter,:) = ...
         [ e_conn(n_el, 3)  e_conn(n_el, 6)  e_conn(n_el, 1) ];
      end


      % A 3-FACE
      if ( tri_adj(n_el,3) < 0 )
        boundary_faces_counter = boundary_faces_counter+1;
        boundary_faces(boundary_faces_counter,:) = ...
         [ e_conn(n_el, 1)  e_conn(n_el, 4)  e_conn(n_el, 2) ];
      end

    end % element loop

    % end quadratic triangle case

  else

    error('tri_mesh_boundary: the routine only accounts for 3 noded\n')
    error('or 6/7 noded triangles\n')

  end

  boundary_faces = boundary_faces(1:boundary_faces_counter,:);

end % function define_boundary_faces
