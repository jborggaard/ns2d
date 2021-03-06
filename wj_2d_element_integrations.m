function [II, JJ, XX, res, entry_counter] = wj_2d_element_integrations(...
                                                                       ...
                                                  x, e_conn,           ...
                                                  ide_u, ide_p,        ...
                                                  u, v,      p,        ...
                                                  uc, vc,   pc,        ...
                                                  material,            ...
                                                  n_pt,                ...
                                                  elem_segment,        ...
                                                  max_part_size,       ...
                                                  tc, t                   )

  %  Integrate the weak residual and Jacobian for the Navier-Stokes equations
  %  over a group of elements.  We assume a Taylor Hood finite element pair.
  
  mu         = material.mu;
  epsilon    = material.epsilon;
  f_function = material.f_function;
  dt         = material.dt;
  theta      = material.theta;
  
  omt        = 1-theta;
  
  
  % use a 7 point rule that handles quintics exactly, for example it can
  % handle terms such as ( u u_x, phi ) = \int (quadratic)*(linear)*(quadratic)
  % a higher order rule should be used for complex forcing terms
  [rr,ss,wt]  = twod_gauss(7);
  one         = ones(size(wt));
  n_equations = max( max(max(ide_u)), max(ide_p) );
  
  % preallocating space
  II        = zeros(15*15*max_part_size,1); % 15 = 2*6 + 1*3
  JJ        = zeros(15*15*max_part_size,1);
  XX        = zeros(15*15*max_part_size,1);
  
  res       = zeros (n_equations, 1);       % weak residual vector contributions
      
  entry_counter = 0;
  for n_el=elem_segment(n_pt)+1:elem_segment(n_pt+1)
    % compute value of each test function and their spatial derivaties
    % at the integration points.
    nodes_local                = e_conn(n_el,:);
    x_local                    = x(nodes_local,:);
    [ ~ , ~  ,phi,phi_x,phi_y] = twod_shape(x_local       ,rr,ss,wt);
    [x_g,wt_g,psi,  ~  ,  ~  ] = twod_shape(x_local(1:3,:),rr,ss,wt);
      
    [fx_g , fy_g ]               = feval(f_function,x_g,t ,mu);
    [fxc_g, fyc_g]               = feval(f_function,x_g,tc,mu);
    % The feval line above must be replaced with the following (or a call to a
    % specific forcing function m-file) to enable Matlab coder to generate a
    % mex file:
    %[fx_g , fy_g ]             = f_function00_2d(x_g,t ,mu);
    %[fxc_g, fyc_g]             = f_function00_2d(x_g,tc,mu);
    
    % compute the value of each variable at the Gauss points
    mu_g      = mu     *one;
    epsilon_g = epsilon*one;
          
    uc_g      = phi  *uc(nodes_local     );
    uc_xg     = phi_x*uc(nodes_local     );
    uc_yg     = phi_y*uc(nodes_local     );
      
    vc_g      = phi  *vc(nodes_local     );
    vc_xg     = phi_x*vc(nodes_local     );
    vc_yg     = phi_y*vc(nodes_local     );
      
    pc_g      = psi  *pc(nodes_local(1:3));
      
    u_g       = phi  *u (nodes_local     );
    u_xg      = phi_x*u (nodes_local     );
    u_yg      = phi_y*u (nodes_local     );
      
    v_g       = phi  *v (nodes_local     );
    v_xg      = phi_x*v (nodes_local     );
    v_yg      = phi_y*v (nodes_local     );
      
    p_g       = psi  *p (nodes_local(1:3));
      
    %---------------------------------------------------------------------------
    %  Assemble the weak form of the equations (element contributions)
    %---------------------------------------------------------------------------
    a11_loc = twod_bilinear(       one, phi  , phi  , wt_g) ...
            + dt*theta*(                                    ...
              twod_bilinear(    2*mu_g, phi_x, phi_x, wt_g) ...  
            + twod_bilinear(      mu_g, phi_y, phi_y, wt_g) ...
            + twod_bilinear(       u_g, phi_x, phi  , wt_g) ...
            + twod_bilinear(       v_g, phi_y, phi  , wt_g) ...
            + twod_bilinear(      u_xg, phi  , phi  , wt_g) ) ;
          
    a12_loc = dt*theta*(                                    ...
              twod_bilinear(      mu_g, phi_x, phi_y, wt_g) ...
            + twod_bilinear(      u_yg, phi  , phi  , wt_g) ) ;
          
    a21_loc = dt*theta*(                                    ...
              twod_bilinear(      mu_g, phi_y, phi_x, wt_g) ...
            + twod_bilinear(      v_xg, phi  , phi  , wt_g) ) ;
          
    a22_loc = twod_bilinear(       one, phi  , phi  , wt_g) ...
            + dt*theta*(                                    ...
              twod_bilinear(      mu_g, phi_x, phi_x, wt_g) ...
            + twod_bilinear(    2*mu_g, phi_y, phi_y, wt_g) ...
            + twod_bilinear(       u_g, phi_x, phi  , wt_g) ...
            + twod_bilinear(       v_g, phi_y, phi  , wt_g) ...
            + twod_bilinear(      v_yg, phi  , phi  , wt_g) ) ;
          
    b1_loc  = twod_bilinear(       one, phi_x, psi  , wt_g)   ;
    b2_loc  = twod_bilinear(       one, phi_y, psi  , wt_g)   ;

    m_loc   = twod_bilinear( epsilon_g, psi  , psi  , wt_g)   ;
        
    f1_loc  =         twod_f_int (                u_g - uc_g , phi  , wt_g)  ...
          + dt*theta*(twod_f_int (     u_g.*u_xg + v_g.*u_yg , phi  , wt_g)  ...
                    - twod_f_int (                       p_g , phi_x, wt_g)  ...
                    + twod_f_int (              2*mu_g.*u_xg , phi_x, wt_g)  ...
                    + twod_f_int (          mu_g.*(u_yg+v_xg), phi_y, wt_g)  ...
                    - twod_f_int (                      fx_g , phi  , wt_g) )...
          + dt*omt  *(twod_f_int ( uc_g.*uc_xg + vc_g.*uc_yg , phi  , wt_g)  ...
                    - twod_f_int (                      pc_g , phi_x, wt_g)  ...
                    + twod_f_int (             2*mu_g.*uc_xg , phi_x, wt_g)  ...
                    + twod_f_int (        mu_g.*(uc_yg+vc_xg), phi_y, wt_g)  ...
                    - twod_f_int (                     fxc_g , phi  , wt_g) )  ;
          
    f2_loc  =         twod_f_int (                v_g - vc_g , phi  , wt_g)  ...
          + dt*theta*(twod_f_int (     u_g.*v_xg + v_g.*v_yg , phi  , wt_g)  ...
                    - twod_f_int (                       p_g , phi_y, wt_g)  ...
                    + twod_f_int (         mu_g.*(u_yg+v_xg) , phi_x, wt_g)  ...
                    + twod_f_int (              2*mu_g.*v_yg , phi_y, wt_g)  ...
                    - twod_f_int (                      fy_g , phi  , wt_g) )...
          + dt*omt  *(twod_f_int ( uc_g.*vc_xg + vc_g.*vc_yg , phi  , wt_g)  ...
                    - twod_f_int (                      pc_g , phi_y, wt_g)  ...
                    + twod_f_int (       mu_g.*(uc_yg+vc_xg) , phi_x, wt_g)  ...
                    + twod_f_int (             2*mu_g.*vc_yg , phi_y, wt_g)  ...
                    - twod_f_int (                     fyc_g , phi  , wt_g) )  ;
          
    f3_loc  =         twod_f_int (               (u_xg + v_yg), psi  , wt_g) ...
            + epsilon*twod_f_int (                        p_g , psi  , wt_g)   ;
        
        
    %---------------------------------------------------------------------------
    %  Assemble contributions into the global system matrix
    %---------------------------------------------------------------------------
    % u-momentum equations
    for n_t=1:6
      n_test = ide_u(nodes_local(n_t),1);
      if (n_test > 0)  % this is an unknown, fill the row
        for n_u=1:6
          n_unku = ide_u(nodes_local(n_u),1);
          n_unkv = ide_u(nodes_local(n_u),2);
          if (n_unku > 0)
            entry_counter = entry_counter + 1;
            II( entry_counter ) = n_test;
            JJ( entry_counter ) = n_unku;
            XX( entry_counter ) = a11_loc(n_t,n_u);
          end
          %
          if (n_unkv > 0)
            entry_counter = entry_counter + 1;
            II( entry_counter ) = n_test;
            JJ( entry_counter ) = n_unkv;
            XX( entry_counter ) = a12_loc(n_t,n_u);
          end
        end
        for n_p=1:3
          n_unkp = ide_p(nodes_local(n_p));

          entry_counter = entry_counter + 1;
          II( entry_counter ) = n_test;
          JJ( entry_counter ) = n_unkp;
          XX( entry_counter ) =-dt*theta*b1_loc(n_p,n_t);
        end
        res(n_test) = res(n_test) + f1_loc(n_t);
      end
    end

      
    % v-momentum equations
    for n_t=1:6
      n_test = ide_u(nodes_local(n_t),2);
      if (n_test > 0)  % this is an unknown, fill the row
        for n_u=1:6
          n_unku = ide_u(nodes_local(n_u),1);
          n_unkv = ide_u(nodes_local(n_u),2);
          if (n_unku > 0)
            entry_counter = entry_counter + 1;
            II( entry_counter ) = n_test;
            JJ( entry_counter ) = n_unku;
            XX( entry_counter ) = a21_loc(n_t,n_u);
          end
          %
          if (n_unkv > 0)
            entry_counter = entry_counter + 1;
            II( entry_counter ) = n_test;
            JJ( entry_counter ) = n_unkv;
            XX( entry_counter ) = a22_loc(n_t,n_u);
          end
        end
        for n_p=1:3
          n_unkp = ide_p(nodes_local(n_p));

          entry_counter = entry_counter + 1;
          II( entry_counter ) = n_test;
          JJ( entry_counter ) = n_unkp;
          XX( entry_counter ) =-dt*theta*b2_loc(n_p,n_t);
        end
        res(n_test) = res(n_test) + f2_loc(n_t);
      end
    end
      
      
    % continuity equations,
    % multiplied by -dt for Jacobian stability, a la Hairer, Lubich, and Roche.
    for n_t=1:3
      n_test = ide_p(nodes_local(n_t));
      for n_u=1:6
        n_unku = ide_u(nodes_local(n_u),1);
        n_unkv = ide_u(nodes_local(n_u),2);
        if (n_unku > 0)
          entry_counter = entry_counter + 1;
          II( entry_counter ) = n_test;
          JJ( entry_counter ) = n_unku;
          XX( entry_counter ) =-dt*b1_loc(n_t,n_u);
        end
        %
        if (n_unkv > 0)
          entry_counter = entry_counter + 1;
          II( entry_counter ) = n_test;
          JJ( entry_counter ) = n_unkv;
          XX( entry_counter ) =-dt*b2_loc(n_t,n_u);
        end
      end
      for n_p=1:3
        n_unkp = ide_p(nodes_local(n_p));

        entry_counter = entry_counter + 1;
        II( entry_counter ) = n_test;
        JJ( entry_counter ) = n_unkp;
        XX( entry_counter ) =-dt*m_loc(n_t,n_p);
      end
      res(n_test) = res(n_test) - dt*f3_loc(n_t);
    end

%       n_uu = find( ide_u(nodes_local,1)>0 );
%       n_vv = find( ide_u(nodes_local,2)>0 );
% 
%       n_unkuu = ide_u(nodes_local(n_uu),1);
%       n_unkvv = ide_u(nodes_local(n_vv),2);
%       n_unkpp = ide_p(nodes_local(1:3));      
%      
%       % u-momentum equations
%       J(n_unkuu, n_unkuu) = J(n_unkuu, n_unkuu) + a11_loc(n_uu,n_uu);
%       J(n_unkuu, n_unkvv) = J(n_unkuu, n_unkvv) + a12_loc(n_uu,n_vv);
%       J(n_unkuu, n_unkpp) = J(n_unkuu, n_unkpp) + b1_loc(1:3,n_uu)';
%       
%       res(n_unkuu) = res(n_unkuu) + f1_loc(n_uu);
% 
%       % v-momentum equations
%       J(n_unkvv, n_unkuu) = J(n_unkvv, n_unkuu) + a21_loc(n_vv,n_uu);
%       J(n_unkvv, n_unkvv) = J(n_unkvv, n_unkvv) + a22_loc(n_vv,n_vv);
%       J(n_unkvv, n_unkpp) = J(n_unkvv, n_unkpp) + b2_loc(1:3,n_vv)';
%       
%       res(n_unkvv) = res(n_unkvv) + f2_loc(n_vv);
% 
%       % continuity equations
%       J(n_unkpp, n_unkuu) = J(n_unkpp, n_unkuu) + b1_loc(1:3,n_uu);
%       J(n_unkpp, n_unkvv) = J(n_unkpp, n_unkvv) + b2_loc(1:3,n_vv);
%       J(n_unkpp, n_unkpp) = J(n_unkpp, n_unkpp) + m_loc;
%       
%       res(n_unkpp) = res(n_unkpp) + f3_loc;
%       

  end % element loop
 

end % function element_integrations


function [fx,fy,q] = f_function00_2d(x,~,~)             %#ok
%% -----------------------------------------------------------------------------
%  Forcing function for the 2D Navier-Stokes equations
%     by convention, the input arguments are: x, t, mu, parameters...
%     the function is called with several quadrature points at once
%     so vector conventions or explicit treatment for each point must be ensured
%-------------------------------------------------------------------------------

  npgaus = size(x,1);
  
  fx = zeros(npgaus,1);
  fy = zeros(npgaus,1);
  q  = zeros(npgaus,1);

end
