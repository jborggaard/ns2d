function [J,res] = weakjac_2d_steady(x, e_conn,    ...
                                     material,     ...
                                     ide_u, ide_p, ...
                                     u, v, p          )
%-------------------------------------------------------------------------------
%  weakjac_2d_steady.m - computes the weak residual of the steady-state
%                        Navier-Stokes equations 
%                        (builds the Newton linearization corresponding to
%                         the discrete Oseen operator when provided u,v,p).
%-------------------------------------------------------------------------------

  f_function  = material.f_function;
  if_mex      = material.if_mex;
  
  n_elements  = size(e_conn,1);

  n_equations = max( max(max(ide_u)), max(ide_p) );
  
  
  %  Integrate finite element matrices
  max_elem_per_partition = 10000;        % maximum number of elements in a mesh
                                         % per partition (for memory efficiency)

  %  Define mesh partitions for integration
  n_part        = floor( n_elements / (max_elem_per_partition+1) ) + 1;
  elem_segment  = floor( linspace(0,n_elements,n_part+1) );

  max_part_size = max( diff( elem_segment ) );
        
  res  = zeros (n_equations, 1); % weak residual vector

  if ( if_mex )
    e_conn = int32(e_conn);
    ide_u  = int32(ide_u);
    ide_p  = int32(ide_p);
    f = functions(f_function);
    if ( ~strcmp(f.function,'f_function0_2d') )  % experimental
      warning('Forcing terms other than zero have not been implemented')
    end
    material.f_function = 'zero';
  end
  
  for n_pt=1:n_part
    if ( if_mex )
      [ II, JJ, XX, res_delta, entry_counter ] =                           ...
                                wj_2d_st_element_integrations_mex(         ...
                                          'wj_2d_st_element_integrations', ...
                                          x, e_conn,                       ...
                                          ide_u, ide_p,                    ...
                                          u, v,      p,                    ...
                                          material,                        ...
                                          int32(n_pt),                     ...
                                          int32(elem_segment),             ...
                                          int32(max_part_size)                );
            
    else
      [ II, JJ, XX, res_delta, entry_counter ] =                           ...
                                     wj_2d_st_element_integrations(        ...
                                          x, e_conn,                       ...
                                          ide_u, ide_p,                    ...
                                          u, v,      p,                    ...
                                          material,                        ...
                                          n_pt,                            ...
                                          elem_segment,                    ...
                                          max_part_size                       );
    end
      
    res = res + res_delta;
    
    if ( n_pt==1 )
      J = sparse( II(1:entry_counter), JJ(1:entry_counter), ...
                  XX(1:entry_counter),                      ...
                  n_equations,         n_equations             );
    else
      J = J + sparse( II(1:entry_counter), JJ(1:entry_counter), ...
                      XX(1:entry_counter),                      ...
                      n_equations,         n_equations             );
    end
    
  end % element partition loop
    
  clear II JJ XX res_delta
   
end % function