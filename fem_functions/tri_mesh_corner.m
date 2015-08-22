function [e_conn, adjacency] = tri_mesh_corner( e_conn, adjacency )
%% -----------------------------------------------------------------------------
%  TRI_MESH_CORNER   Alters meshes to rearrange elements w/ two boundary faces.
%
%    So-called "meshing into corners" is important in simulating incompressible 
%    fluid flows.
%
%  Version: 1.2
%
%  Usage:    [e_conn,adjacency] = tri_mesh_corner(e_conn,adjacency);
%
%
%    where 
%           e_conn     is the element connectivity (must have 3 or 6 columns)
%           adjacency  is the element adjacency (optional, if not provided, it
%                      will be computed and is available for return)
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    23 May 2013
%
%  Authors:
%
%    Jeff Borggaard and John Burkardt
%% -----------------------------------------------------------------------------
  verbose = 0;

  [tri_num, tri_order] = size(e_conn);
  
  if ( tri_order ~= 3 && tri_order ~= 6 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'TRI_MESH_CORNER - Fatal error!\n' );
    fprintf ( 1, '  Data is not for a 3-node or 6-node triangulation.\n' );
    error ( 'TRI_MESH_CORNER - Fatal error!' );
  end

%
%  If needed, create the triangle neighbor array.
%
  if ( nargin<3 )
    adjacency = tri_mesh_neighbors ( e_conn );
  end

%
%  Examine the triangle neighbor array.
%
  num_bdy_faces = sum( adjacency<0, 2 );
  
  OFFSET = 1;
  negative_total(OFFSET+0:OFFSET+3) = 0;

  for triangle = 1 : tri_num

    negative_total(OFFSET+num_bdy_faces(triangle)) = ...
        negative_total(OFFSET+num_bdy_faces(triangle)) + 1;

  end

  if ( verbose )
    fprintf ( 1, '\n' );
    fprintf ( 1, '  Number of boundary sides     Number of triangles\n' );
    fprintf ( 1, '\n' );
    for i = 0 : 3
      fprintf ( 1, '            %8d            %8d\n', i, negative_total(OFFSET+i) );
    end
  end
%
%  Try to patch problems.
%
  if ( 0 < negative_total(OFFSET+3) )

    fprintf ( 1, '\n' );
    fprintf ( 1, 'TRI_MESH_CORNER - Fatal error!\n' );
    fprintf ( 1, '  There is at least one triangle with all sides\n' );
    fprintf ( 1, '  on the boundary.\n' );
    return

  elseif ( 0 == negative_total(OFFSET+2) )

    fprintf ( 1,  '\n' );
    fprintf ( 1,  'TRI_MESH_CORNER:\n' );
    fprintf ( 1,  '  No corner triangles were found.\n' );
    fprintf ( 1,  '  No corrections need to be made.\n' );
    return

  else
    %
    %  Consider those triangles with exactly two boundary sides.
    %
    for triangle1 = 1 : tri_num 

      if ( num_bdy_faces(triangle1) == 2 )

        triangle2 = -1;

        for neighbor = 1 : 3
          if ( 0 < adjacency(triangle1,neighbor) )
            triangle2 = adjacency(triangle1,neighbor);
            t1_to_t2 = neighbor;
          end
        end

        if ( verbose )
          fprintf ( 1, '  Adjusting triangle %d using triangle %d\n', ...
            triangle1, triangle2 );
        end
        
        t2_to_t1 = -1;
        for neighbor = 1 : 3
          if ( adjacency(triangle2,neighbor) == triangle1 )
            t2_to_t1 = neighbor;
          end
        end

        nodes1 = e_conn(triangle1,:);
        nodes2 = e_conn(triangle2,:);

        if ( tri_order == 3 )
          switch (t1_to_t2)
            case 1 
            case 2, nodes1 = nodes1( [2 3 1] );
            case 3, nodes1 = nodes1( [3 1 2] );
          end
          
          switch (t2_to_t1)
            case 1
            case 2, nodes2 = nodes2( [2 3 1] );
            case 3, nodes2 = nodes2( [3 1 2] );
          end
          
          e_conn( triangle1,: ) = [ nodes1(1) nodes2(1) nodes1(3) ];
          e_conn( triangle2,: ) = [ nodes1(1) nodes1(2) nodes2(1) ];
          
        else
          switch (t1_to_t2)
            case 1 
            case 2, nodes1 = nodes1( [2 3 1 5 6 4] );
            case 3, nodes1 = nodes1( [3 1 2 6 4 5] );
          end
          
          switch (t2_to_t1)
            case 1
            case 2, nodes2 = nodes2( [2 3 1 5 6 4] );
            case 3, nodes2 = nodes2( [3 1 2 6 4 5] );
          end
          
          e_conn( triangle1,: ) = [ nodes1(1) nodes2(1) nodes1(3) ...
                                    nodes1(5) nodes2(4) nodes1(6)    ];
          e_conn( triangle2,: ) = [ nodes1(1) nodes1(2) nodes2(1) ...
                                    nodes1(4) nodes2(6) nodes2(5)    ];
        end
        
        %
        %  Update the adjacency array.
        %
        switch (t2_to_t1)
          case 1, tri2 = adjacency(triangle2,2);
                  tri3 = adjacency(triangle2,3);
          case 2, tri2 = adjacency(triangle2,3);
                  tri3 = adjacency(triangle2,1);
          case 3, tri2 = adjacency(triangle2,1);
                  tri3 = adjacency(triangle2,2);
        end
        
        adjacency(triangle1,:) = [ tri3, -1, triangle2 ];
        adjacency(triangle2,:) = [ tri2, triangle1, -1 ];

        for i=1:3
          if ( adjacency( tri3, i ) == triangle2 )
            adjacency(tri3,i) = triangle1;
          end
        end
        
      end

    end

  end
  
end