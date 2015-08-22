function [neighbors,v_conn] = tri_mesh_neighbors ( e_conn )
%%--------------------------------------------------------------------------78--
%
%  MAIN is the main program for TRI_MESH_NEIGHBORS.
%
%  Discussion:
%
%    TRI_MESH_NEIGHBORS determines the neighboring elements for each triangular
%    element and for each vertex in a triangular mesh.  This is based on the
%    function TRIANGULARION_TRIANGLE_NEIGHBORS written by John Burkardt.
%
%  Usage:
%
%    [neighbors,v_conn] = tri_mesh_neighbors ( e_conn );
%
%    where:
%
%    e_conn     is the element connectivity.  Standard orientation is assumed
%               so the first three nodes in the connectivity are the vertices.
%               [ n_elements, n_local_nodes ]
%
%    neighbors  is a [ n_elements, 3 ] array which contains the adjacent 
%               element to each triangle edge.  A value of -1 is given if that
%               edge has no adjacent element (e.g., on the boundary)
%
%                                3
%                               | \
%                        edge 2 |  \ edge 1
%                               | e \
%                               1----2
%                               edge 3
%
%
%    v_conn     is a structure of n_node arrays that list elements that are
%               connected to each node.  v_conn(1).e_list contains a list of the
%               elements that include node 1, etc.  ( optional output )
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    09 June 2012
%
%  Authors:
%
%    Jeff Borggaard and John Burkardt
%
%%--------------------------------------------------------------------------78--

  verbose = 0; 
  if ( verbose )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'TRI_MESH_NEIGHBORS\n' );
  end
  
  if ( nargin < 1 )
    error ('TRI_MESH_NEIGHBORS: no connectivity entered' );
  end

  [ n_elements, nel_dof ] = size( e_conn );
  
%
%  Detect and correct 0-based indexing.
%
  min_n = min( min( e_conn ) );
  
  if ( min_n==0 )
    fprintf('TRI_MESH_NEIGHBORS: element connectivity may have 0-based ')
    fprintf('indexing\n will attempt to adjust and proceed\n')
    e_conn = e_conn + 1;
  end
  
  n_nodes = max( max( e_conn ) );
  
%
%  Create the triangle neighbor array.
%
%
%  Step 1.
%  From the list of nodes for triangle T, of the form: (I,J,K)
%  construct the three neighbor relations:
%
%    (I,J,3,T) or (J,I,3,T),
%    (J,K,1,T) or (K,J,1,T),
%    (K,I,2,T) or (I,K,2,T)
%
%  where we choose (I,J,3,T) if I < J, or else (J,I,3,T)
%
  if ( verbose )
    fprintf('TRI_MESH_NEIGHBORS: defining sorted faces\n');
  end
  
  faces = zeros(4,3*n_elements);
  for nel = 1 : n_elements

    i = e_conn(nel,1);
    j = e_conn(nel,2);
    k = e_conn(nel,3);

    if ( i < j )
      faces(1:4,1+3*(nel-1)) = [ i, j, 3, nel ]';
    else
      faces(1:4,1+3*(nel-1)) = [ j, i, 3, nel ]';
    end

    if ( j < k )
      faces(1:4,2+3*(nel-1)) = [ j, k, 1, nel ]';
    else
      faces(1:4,2+3*(nel-1)) = [ k, j, 1, nel ]';
    end

    if ( k < i )
      faces(1:4,3+3*(nel-1)) = [ k, i, 2, nel ]';
    else
      faces(1:4,3+3*(nel-1)) = [ i, k, 2, nel ]';
    end

  end
%
%  Step 2. Perform an ascending dictionary sort on the neighbor relations.
%  We only intend to sort on rows 1:2; the routine we call here
%  sorts on rows 1 through 4 but that won't hurt us.
%
%  What we need is to find cases where two triangles share a face.
%  By sorting the columns of the FACES array, we will put shared faces
%  next to each other.
%
%  this scrambles the element information.  
%   faces = sort(faces,1);
   
   if ( verbose )
     fprintf('TRI_MESH_NEIGHBORS: grouping faces\n');
   end
   
   [~,index2] = sort(faces(2,:));
   faces = faces(:,index2);
   
   [~,index1] = sort(faces(1,:));
   faces = faces(:,index1);
%
%  Step 3. Neighboring triangles show up as consecutive columns with
%  identical first two entries.  Whenever you spot this happening,
%  make the appropriate entries in TRIANGLE_NEIGHBOR.
%
  neighbors= -1*ones(n_elements,3);

  if ( verbose )
    fprintf('TRI_MESH_NEIGHBORS: finding matching faces\n');
  end
  
  face = 1;
  while ( 1 )

    if ( 3 * n_elements <= face )
      break
    end

    if ( faces(1:2,face) == faces(1:2,face+1) )
      face1 = faces(3,face);
      elem1 = faces(4,face);
      face2 = faces(3,face+1);
      elem2 = faces(4,face+1);
      neighbors(elem1,face1) = elem2;
      neighbors(elem2,face2) = elem1;
      face = face + 2;
    else % there is no neighbor, move on
      face = face + 1;
    end

  end
    
  
  if ( nargout==2 )
    for n=1:n_nodes
      v_conn(n).e_list = [];
    end
    for n=1:n_elements  % for each node, add the current element to e_list
      for m=1:nel_dof
      v_conn(e_conn(n,m)).e_list = [ v_conn(e_conn(n,m)).e_list, n ];
    end
  end
end

end % function tri_mesh_neighbors