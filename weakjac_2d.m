function [J,res] = weakjac_2d(x, e_conn,    ...
                              material,     ...
                              ide_u, ide_p, ...
                              u , v , p ,   ...
                              uc, vc, pc,   ...
                              tc, t            )
%% -----------------------------------------------------------------------------
%  weakjac_2d.m - computes the weak residual of the 2d Navier-Stokes
%                 equations.
%
%  Variables:
%                 x           nodal coordinates
%                 e_conn      element connectivity
%                 material    structure containing flow and solution parameters
%                             mu - viscosity (1/Re)
%                             epsilon - penalty parameter
%                             f_function - pointer to forcing function
%                             if_mex     - flag for mex element integration
%                             dt         - time step
%                             theta      - numerical integration parameter
%
%                 ide_u       velocity equation numbers for each node
%                 ide_p       pressure equation numbers for each node
%                 u,v,p       Newton iterate of solution at next timestep
%                 uc,vc,pc    solution at current timestep
%                 tc,t        current and next time
%%------------------------------------------------------------------------------

  f_function  = material.f_function;
  if_mex      = material.if_mex;
  
  n_elements  = size(e_conn,1);

  n_equations = max( max(max(ide_u)), max(ide_p) );
  
  
  %  Integrate finite element matrices
  max_elem_per_partition = 10000;  % maximum number of elements in a mesh
                                   % per partition (for memory efficiency)

  %  Define mesh partitions for integration
  n_part        = floor( n_elements / (max_elem_per_partition+1) ) + 1;
  elem_segment  = floor( linspace(0,n_elements,n_part+1) );
%  fprintf('  Partitioned elements into segments:\n')
%  disp(elem_segment)
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
      [ II, JJ, XX, res_delta, entry_counter ] = ...
      wj_2d_element_integrations_mex(            ...
                'wj_2d_element_integrations',    ...
                x, e_conn,                       ...
                ide_u, ide_p,                    ...
                u, v,      p,                    ...
                uc, vc,   pc,                    ...
                material,                        ...
                int32(n_pt),                     ...
                int32(elem_segment),             ...
                int32(max_part_size),            ...
                tc, t                               );
            
    else
      [ II, JJ, XX, res_delta, entry_counter ] = ...
      wj_2d_element_integrations(                ...
                x, e_conn,                       ...
                ide_u, ide_p,                    ...
                u, v,      p,                    ...
                uc, vc,   pc,                    ...
                material,                        ...
                n_pt,                            ...
                elem_segment,                    ...
                max_part_size,                   ...
                tc, t                               );
    end
      
%    fprintf('  finished integrating partition segment %d\n',n_pt);
    
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
    
  end % part
    
  clear II JJ XX res_delta
   
end % function
