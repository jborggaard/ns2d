function [d1_p, d2_p, e_error, node] = twod_projectd(x,e_conn,u,node)
%-------------------------------------------------------------------------------
%  twod_projectd.m - routine to project the derivative of scalar fields onto 
%                    the continuous finite element space (ZZ projection).
%
%  Copyright (c) 2001, Jeff Borggaard, Virginia Tech
%  Version: 2.0
%
%  Usage:   [dudx1, dudx2, e_error, node] = twod_projectd(x, e_conn, u, node)
%
%                 Multiple scalar fields can be treated at once 
%                     u = [u1  u2  ... ]
% 
%                  in which case  
%                     dudx1 = [d(s1)/d(x1)  d(s2)/d(x1)  ... ]
%                     dudx2 = [d(s1)/d(x1)  d(s2)/d(x1)  ... ]
%
%  Output Variables:     
%                 dudx1   - Projection of the x-derivative evaluated at nodes
%                 dudx2   - Projection of the y-derivative evaluated at nodes
%
%                 e_error - element error (H1-seminorm for each field)
%                 node    - a structure that contains an element list for
%                           each node
%  Input Variables:    
%                 x       - nodal coordinates
%                 e_conn  - element connectivity
%                 u       - nodal values of scalar quantity
%
%                 node    - (optional)
%                           the node structure can be reused for multiple
%                            projections involving the same mesh.
%% -----------------------------------------------------------------------------

  n_nodes    = size(x,1);
  n_elements = size(e_conn,1);

  n_fields   = size(u,2);

  %-------------------------------------------------------------------------------
  %  For every vertex node, construct a list of elements that share it
  %-------------------------------------------------------------------------------
  node(n_nodes).e_list = [];

  for n=1:n_elements  % local nodes 1,2, and 3 are vertices
    node(e_conn(n,1)).e_list = [ node(e_conn(n,1)).e_list, n ];
    node(e_conn(n,2)).e_list = [ node(e_conn(n,2)).e_list, n ];
    node(e_conn(n,3)).e_list = [ node(e_conn(n,3)).e_list, n ];
  end

  %-------------------------------------------------------------------------------
  %  Compute the least-squares projection over each patch,
  %     least-squares of gradients are performed at element quadrature points
  %-------------------------------------------------------------------------------
  n_gauss = 7;
  [r,s,w] = twod_gauss(n_gauss);

  d1_p = realmax*ones(n_nodes,n_fields);
  d2_p = realmax*ones(n_nodes,n_fields);  % initialize to ridiculus values

  for n=1:n_nodes
    if ( ~isempty(node(n).e_list) ) % this is a vertex node
      nel = length(node(n).e_list);
      P   = zeros(nel*n_gauss,6); 
      d1  = zeros(nel*n_gauss,n_fields);
      d2  = d1;
    
      for j=1:length(node(n).e_list) % loop over elements in the patch
        el = node(n).e_list(j);
        nodes_local = e_conn(el,:);
        x_local = x(nodes_local,:);
        u_local = u(nodes_local,:); 

        [x_g, ~ , ~ ,p_x, p_y] = twod_shape(x_local, r, s, w);

        d_x = p_x*u_local;   % calculate the finite element derivatives
        d_y = p_y*u_local;   % on this element

        idx = 1+(j-1)*n_gauss:j*n_gauss;
        P(idx,:) = ...
        [ ones(n_gauss,1) x_g(:,1) x_g(:,2) x_g(:,1).^2 x_g(:,1).*x_g(:,2) x_g(:,2).^2 ];
        d1(idx,:) = d_x;
        d2(idx,:) = d_y;
      
%       for i=1:n_gauss
%         P  = [ P ; ...
%                1 x_g(i,1) x_g(i,2) x_g(i,1)^2 x_g(i,1)*x_g(i,2) x_g(i,2)^2 ];
%         d1 = [ d1; d_x(i) ];
%         d2 = [ d2; d_y(i) ];
%       end
      end % finished filling least-squares system

      % calculate the polynomial coefficients for each field
      for i=1:n_fields
        a_x = P\d1(:,i);
        a_y = P\d2(:,i);

        %-------------------------------------------------------------------------
        %  Compute projected derivatives at vertex node
        %-------------------------------------------------------------------------
        xp  = x(n,1);
        yp  = x(n,2);
        d1_p(n,i) = a_x(1)      + a_x(2)*xp    + a_x(3)*yp + ...
                    a_x(4)*xp^2 + a_x(5)*xp*yp + a_x(6)*yp^2;
        d2_p(n,i) = a_y(1)      + a_y(2)*xp    + a_y(3)*yp + ...
                    a_y(4)*xp^2 + a_y(5)*xp*yp + a_y(6)*yp^2;

        for j=1:length(node(n).e_list)
          el = node(n).e_list(j);
          nodes_local = e_conn(el,:);
          x_local = x(nodes_local,:);

          %-----------------------------------------------------------------------
          %  Compute contribution to edge nodes.  These are either on the
          %  boundary or are averaged over two element faces.
          %-----------------------------------------------------------------------
          if ( n==nodes_local(1) )
            x4 = x_local(4,1);  y4 = x_local(4,2);
            x6 = x_local(6,1);  y6 = x_local(6,2);

            if ( d1_p(nodes_local(4),i) > realmax/2 ) % contribution not rec'd
              d1_p(nodes_local(4),i) = ...
                  a_x(1)      + a_x(2)*x4    + a_x(3)*y4 + ...
                  a_x(4)*x4^2 + a_x(5)*x4*y4 + a_x(6)*y4^2;
              d2_p(nodes_local(4),i) = ...
                  a_y(1)      + a_y(2)*x4    + a_y(3)*y4 + ...
                  a_y(4)*x4^2 + a_y(5)*x4*y4 + a_y(6)*y4^2;
            else
              d1_p(nodes_local(4),i) = .5*( d1_p(nodes_local(4),i) + ...
                  a_x(1)      + a_x(2)*x4    + a_x(3)*y4 + ...
                  a_x(4)*x4^2 + a_x(5)*x4*y4 + a_x(6)*y4^2 );
              d2_p(nodes_local(4),i) = .5*( d2_p(nodes_local(4),i) + ...
                  a_y(1)      + a_y(2)*x4    + a_y(3)*y4 + ...
                  a_y(4)*x4^2 + a_y(5)*x4*y4 + a_y(6)*y4^2 );
            end
            if ( d1_p(nodes_local(6),i) > realmax/2 )
              d1_p(nodes_local(6),i) = ...
                  a_x(1)      + a_x(2)*x6    + a_x(3)*y6 + ...
                  a_x(4)*x6^2 + a_x(5)*x6*y6 + a_x(6)*y6^2;
              d2_p(nodes_local(6),i) = ...
                  a_y(1)      + a_y(2)*x6    + a_y(3)*y6 + ...
                  a_y(4)*x6^2 + a_y(5)*x6*y6 + a_y(6)*y6^2;
            else
              d1_p(nodes_local(6),i) = .5*( d1_p(nodes_local(6),i) + ...
                  a_x(1)      + a_x(2)*x6    + a_x(3)*y6 + ...
                  a_x(4)*x6^2 + a_x(5)*x6*y6 + a_x(6)*y6^2 );
              d2_p(nodes_local(6),i) = .5*( d2_p(nodes_local(6),i) + ...
                  a_y(1)      + a_y(2)*x6    + a_y(3)*y6 + ...
                  a_y(4)*x6^2 + a_y(5)*x6*y6 + a_y(6)*y6^2 );
            end

          elseif ( n==nodes_local(2) ) 
            x4 = x_local(4,1);  y4 = x_local(4,2);
            x5 = x_local(5,1);  y5 = x_local(5,2);
  
            if ( d1_p(nodes_local(4),i) > realmax/2 )
              d1_p(nodes_local(4),i) = ...
                  a_x(1)      + a_x(2)*x4    + a_x(3)*y4 + ...
                  a_x(4)*x4^2 + a_x(5)*x4*y4 + a_x(6)*y4^2;
              d2_p(nodes_local(4),i) = ...
                  a_y(1)      + a_y(2)*x4    + a_y(3)*y4 + ...
                  a_y(4)*x4^2 + a_y(5)*x4*y4 + a_y(6)*y4^2;
            else
              d1_p(nodes_local(4),i) = .5*( d1_p(nodes_local(4),i) + ...
                  a_x(1)      + a_x(2)*x4    + a_x(3)*y4 + ...
                  a_x(4)*x4^2 + a_x(5)*x4*y4 + a_x(6)*y4^2 );
              d2_p(nodes_local(4),i) = .5*( d2_p(nodes_local(4),i) + ...
                  a_y(1)      + a_y(2)*x4    + a_y(3)*y4 + ...
                  a_y(4)*x4^2 + a_y(5)*x4*y4 + a_y(6)*y4^2 );
            end
            if ( d1_p(nodes_local(5),i) > realmax/2 )
              d1_p(nodes_local(5),i) = ...
                  a_x(1)      + a_x(2)*x5    + a_x(3)*y5 + ...
                  a_x(4)*x5^2 + a_x(5)*x5*y5 + a_x(6)*y5^2;
              d2_p(nodes_local(5),i) = ...
                  a_y(1)      + a_y(2)*x5    + a_y(3)*y5 + ...
                  a_y(4)*x5^2 + a_y(5)*x5*y5 + a_y(6)*y5^2;
            else
              d1_p(nodes_local(5),i) = .5*( d1_p(nodes_local(5),i) + ...
                  a_x(1)      + a_x(2)*x5    + a_x(3)*y5 + ...
                  a_x(4)*x5^2 + a_x(5)*x5*y5 + a_x(6)*y5^2 );
              d2_p(nodes_local(5),i) = .5*( d2_p(nodes_local(5),i) + ...
                  a_y(1)      + a_y(2)*x5    + a_y(3)*y5 + ...
                  a_y(4)*x5^2 + a_y(5)*x5*y5 + a_y(6)*y5^2 );
            end

          elseif ( n==nodes_local(3) ) 
            x5 = x_local(5,1);  y5 = x_local(5,2);
            x6 = x_local(6,1);  y6 = x_local(6,2);

            if ( d1_p(nodes_local(5),i) > realmax/2 )
              d1_p(nodes_local(5),i) = ...
                  a_x(1)      + a_x(2)*x5    + a_x(3)*y5 + ...
                  a_x(4)*x5^2 + a_x(5)*x5*y5 + a_x(6)*y5^2;
              d2_p(nodes_local(5),i) = ...
                  a_y(1)      + a_y(2)*x5    + a_y(3)*y5 + ...
                  a_y(4)*x5^2 + a_y(5)*x5*y5 + a_y(6)*y5^2;
            else
              d1_p(nodes_local(5),i) = .5*( d1_p(nodes_local(5),i) + ...
                  a_x(1)      + a_x(2)*x5    + a_x(3)*y5 + ...
                  a_x(4)*x5^2 + a_x(5)*x5*y5 + a_x(6)*y5^2 );
              d2_p(nodes_local(5),i) = .5*( d2_p(nodes_local(5),i) + ...
                  a_y(1)      + a_y(2)*x5    + a_y(3)*y5 + ...
                  a_y(4)*x5^2 + a_y(5)*x5*y5 + a_y(6)*y5^2 );
            end
            if ( d1_p(nodes_local(6),i) > realmax/2 )
              d1_p(nodes_local(6),i) = ...
                  a_x(1)      + a_x(2)*x6    + a_x(3)*y6 + ...
                  a_x(4)*x6^2 + a_x(5)*x6*y6 + a_x(6)*y6^2;
              d2_p(nodes_local(6),i) = ...
                  a_y(1)      + a_y(2)*x6    + a_y(3)*y6 + ...
                  a_y(4)*x6^2 + a_y(5)*x6*y6 + a_y(6)*y6^2;
            else
              d1_p(nodes_local(6),i) = .5*( d1_p(nodes_local(6),i) + ...
                  a_x(1)      + a_x(2)*x6    + a_x(3)*y6 + ...
                  a_x(4)*x6^2 + a_x(5)*x6*y6 + a_x(6)*y6^2 );
              d2_p(nodes_local(6),i) = .5*( d2_p(nodes_local(6),i) + ...
                  a_y(1)      + a_y(2)*x6    + a_y(3)*y6 + ...
                  a_y(4)*x6^2 + a_y(5)*x6*y6 + a_y(6)*y6^2 );
            end
          else
            fprintf('Problem identifying node at a vertex')
          end

        end
      end
    end
  end
%  Compare d1_p with the true gradient
%d_true=-2*x(:,1).*(1-x(:,2).^2); % elliptic_2da example
%disp([d1_p d_true])

%  Compare d2_p with the true gradient
%d_true=-2*x(:,2).*(1-x(:,1).^2);
%disp([d2_p d_true])

  if ( nargout>2 )
    %-----------------------------------------------------------------------------
    %  Calculate the H1-seminorm of the error on each element
    %-----------------------------------------------------------------------------
    e_error = zeros(n_elements,1);

    n_gauss = 7;
    [r,s,w]   = twod_gauss(n_gauss);

    for n_el=1:n_elements
      nodes_local       = e_conn(n_el,:);
      x_local           = x(nodes_local,:);
      u_local           = u(nodes_local,:);
      d1_local          = d1_p(nodes_local,:);
      d2_local          = d2_p(nodes_local,:);

      [~,w_g,phi,p_x,p_y] = twod_shape(x_local,r,s,w);

      % calculate the finite element derivatives and their projections at x_g
      dx_fe           = p_x*u_local;
      dy_fe           = p_y*u_local;

      dx_p            = phi*d1_local;
      dy_p            = phi*d2_local;

  % % exact derivatives as a test
  % dx_p = -2*x_g(:,1).*(1-x_g(:,2).^2);
  % dy_p = -2*x_g(:,2).*(1-x_g(:,1).^2);
      for i=1:n_fields
        sq_error     = (dx_fe(:,i)-dx_p(:,i)).*(dx_fe(:,i)-dx_p(:,i)) + (dy_fe(:,i)-dy_p(:,i)).*(dy_fe(:,i)-dy_p(:,i));
        e_error(n_el,i) = sq_error'*w_g;
  
      end

    e_error = sqrt(e_error);

  end

end
 % min_e = min(e_error);  max_e = max(e_error);  range_e = max_e-min_e;

%   if (nargin==4)
%     figure(fig_num)
%     for el=1:n_elements
%       n1 = e_conn(el,1);  n2 = e_conn(el,2);  n3 = e_conn(el,3);
% 
%       x1 = x(n1,1); y1 = x(n1,2); e1 = (e_error(el)-min_e)/range_e;
%       x2 = x(n2,1); y2 = x(n2,2); e2 = (e_error(el)-min_e)/range_e;
%       x3 = x(n3,1); y3 = x(n3,2); e3 = (e_error(el)-min_e)/range_e;
% 
%       patch([x1 x2 x3]',[y1 y2 y3]',[e1 e2 e3]')
%     end
% 
%     axis equal;
%   end


