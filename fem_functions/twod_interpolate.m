function [f,element_list] = twod_interpolate(      ...
                   xy_node,                                        ...
                   e_conn,                                         ...
                   neighbors,                                      ...
                   f_node,                                         ...
                   xy_points)
%%
%  TWOD_INTERPOLATE  A MATLAB function to interpolate a finite element
%  solution at a collection of points.  This function assumes that the
%  triangular mesh is a Delauney triangulation and will retry searches if
%  interior holes are encountered.
%
%  Useage:
%
%  Variables:
%        xy_node          - finite element nodal coordinates
%        e_conn           - finite element connectivity
%        neighbors        - element adjacency matrix
%        f_node           - values of given functions at the nodes
%
%        xy_points        - interpolation points
%
%        f                - interpolated values of f_node at xy_points
%        element_list     - a list of elements that contain xy_points
%%
%  save('cylinder.node','-ascii',xy)
%
  [n_points,dim] = size(xy_points);
  [n_elem,npe]   = size(e_conn);

  element_list    = zeros(n_points,1);
  iso_coordinates = zeros(n_points,npe);

  for i=1:n_points
    [element_list(i), bary1, bary2] = mesh_search(xy_node,         ...
                                    e_conn(:,1:3),                 ...
                                    neighbors,                     ...
                                    xy_points(i,:));
    
    nodes_local = e_conn( element_list(i),: );
    [x_g,w_g,phi,p_x,p_y] = twod_shape(xy_node(nodes_local,:),     ...
                                       bary1,                      ...
                                       bary2,                      ...
                                       [1]);

    f(i,:) = phi*f_node(nodes_local,:);
  end
   
end

function [ element, iso1, iso2 ] = mesh_search( xy, e_conn, neighbors, point )
%%
%  MESH_SEARCH  A MATLAB function that locates which element contains a
%  given point.  The function only provides one such element and makes an
%  explicit assumption that the triangulation is Delauney.
%%

  [n_elem,npe] = size(e_conn);
%  figure
%  hold on

  element = 1;   % initial guess of the element containing the point
  element_old = 1;
  count = 0;
  found = 0;
  while (~found)
    x_local = xy( e_conn(element,1:3), : );
    
%    plot( [ x_local(1,1) x_local(2,1) ],[ x_local(1,2) x_local(2,2)],...
%          [ x_local(2,1) x_local(3,1) ],[ x_local(2,2) x_local(3,2)],...
%          [ x_local(3,1) x_local(1,1) ],[ x_local(3,2) x_local(1,2)] )
      
    delta = ( x_local(3,1)-x_local(2,1) )*( x_local(1,2)-x_local(2,2) ) ...
          - ( x_local(3,2)-x_local(2,2) )*( x_local(1,1)-x_local(2,1) );
      
    iso1  = ( ( x_local(2,1)-point(1) )*( x_local(3,2)-point(2) ) ...
            - ( x_local(3,1)-point(1) )*( x_local(2,2)-point(2) ) ) / delta;
  
    iso2  = ( ( x_local(3,1)-point(1) )*( x_local(1,2)-point(2) ) ...
            - ( x_local(1,1)-point(1) )*( x_local(3,2)-point(2) ) ) / delta;
  
    iso3  = ( ( x_local(1,1)-point(1) )*( x_local(2,2)-point(2) ) ...
            - ( x_local(2,1)-point(1) )*( x_local(1,2)-point(2) ) ) / delta;
    
    if ( iso1>0 & iso2>0 & iso3>0 )
      found = 1;
    
    else
      count = count + 1;
      [tmp,index] = sort( [iso1 iso2 iso3],'ascend' );
      
      if ( neighbors(element,index(1)) > 0 )
        element_old = element;
        element = neighbors(element,index(1));
        choice = 1;
      elseif ( neighbors(element,index(2)) > 0 )
        element_old = element;
        element = neighbors(element,index(2));
        choice = 2;
      elseif ( neighbors(element,index(3)) > 0 )
        element_old = element;
        element = neighbors(element,index(3));
        choice = 3;
      else
          iso1
          iso2
          iso3
          neighbors(element,:)
          index
          stop
      end
%       if ( index(choice)==1 )
%         plot( (x_local(2,1)+x_local(3,1))/2, (x_local(2,2)+x_local(3,2))/2, '*')
%       elseif ( index(choice)==2 )
%         plot( (x_local(1,1)+x_local(3,1))/2, (x_local(1,2)+x_local(3,2))/2, '*')
%       else
%         plot( (x_local(1,1)+x_local(2,1))/2, (x_local(1,2)+x_local(2,2))/2, '*')
%       end
       if ( count > n_elem )   % if we are somehow in a loop, get out...
         element = n_elem;
         count   = 0;
       end
    end
    
  end
  
%  count
end