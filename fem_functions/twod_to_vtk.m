function [] = twod_to_vtk(filename,x,e_conn,scalars,vectors,labels)
%-------------------------------------------------------------------------------
%  twod_to_vtk.m - Writes out a 2D finite element solution in the 
%                  VTK format (visualization toolkit).  This format
%                  can be used by, e.g. ParaView.
%
%  Copyright (c) 2010, Jeff Borggaard, Virginia Tech
%  Version: 1.1
%
%  Usage:    [] = twod_to_vtk(filename,x,e_conn,scalars,vectors,labels)
%
%                 twod_to_vtk(filename,x,e_conn,p,[],{'pressure'})
%
%                 twod_to_vtk(filename,x,e_conn,[],[u v],{'velocity'})
%
%                 twod_to_vtk(filename,x,e_conn,[p T],[u v],{'p','tmp','vel'})
%
%  Variables:     filename
%                        Name of the output file
%                 x
%                        Nodal coordinates
%                 e_conn
%                        Element connectivity
%                 scalars
%                        Columns of scalar variables to write to vtk file
%                        (can be [])
%                 vectors
%                        Must be ( n_nodes x (2*n_vector sets) )
%                 labels
%                        A cell array with names of scalar and vector variables
%-------------------------------------------------------------------------------

  % Test input files and extract sizes
  %-----------------------------------------------------------------------------
  if ( isempty(filename) )
    error('filename must be specified');
  end
  
  out_vtk = fopen(filename,'w');

  [n_nodes   , n_dim      ] = size(x);
  [n_elements, nel_dof    ] = size(e_conn);
  [sn        , n_scalars  ] = size(scalars);
  [vn        , n_vectors  ] = size(vectors);
 
  if ( n_dim ~= 2 )
    warning('expecting two coordinates, continuing anyway...')
  end
  
  if ( n_scalars > 0 && n_nodes ~= sn )
    error('incompatible dimensions between nodes and scalars\n')
  end
  
  if ( n_vectors > 0 && n_nodes ~= vn )
    error('incompatible dimensions between nodes and vectors\n')
  end

  if ( mod(n_vectors,2) )
    error('vectors must be two dimensional\n')
  end
  
  if ( nargin<6 || length(labels)<n_scalars+n_vectors/2 )
    for i=1:n_scalars
      labels{i} = ['scalar',int2str(i)];
    end
    for i=1:n_vectors/2
      labels{n_scalars+i} = ['vector',int2str(i)];
    end
  end
 
  if ( n_scalars>0 )
    scalars(abs(scalars)<eps)=0;
  end
  
  if ( n_vectors>0 )
    vectors(abs(vectors)<eps)=0;
  end
  
  % Preamble
  %-----------------------------------------------------------------------------
  fprintf(out_vtk,'# vtk DataFile Version 2.0\n');
  fprintf(out_vtk,'Solution: %s, Timestep %d\n','title',0);
  fprintf(out_vtk,'ASCII\n\n');
  fprintf(out_vtk,'DATASET UNSTRUCTURED_GRID\n');
  
  
  % Write out coordinate values
  %-----------------------------------------------------------------------------
  fprintf(out_vtk,'POINTS %8d float\n',n_nodes);
  for i=1:n_nodes
    fprintf(out_vtk,'%12.4e %12.4e %12.4e\n',x(i,1),x(i,2),0);
  end
  
  % Write out element connectivity
  %-----------------------------------------------------------------------------
  if ( nel_dof==3 )                               % Test for linear triangles
    fprintf(out_vtk,'CELLS %7d %7d\n',n_elements,4*n_elements);

    for i=1:n_elements
      fprintf(out_vtk,'3  %6i  %6i  %6i\n',e_conn(i,:)-1);
    end

    fprintf(out_vtk,'\nCELL_TYPES %d\n',n_elements);
    for i=1:n_elements
      fprintf(out_vtk,'5\n');
    end
    
  elseif ( nel_dof==6 )                           % Test for quadratic triangles
    fprintf(out_vtk,'CELLS %7d %7d\n',n_elements,7*n_elements);

    for i=1:n_elements
      fprintf(out_vtk,'6  %6i  %6i  %6i  %6i  %6i  %6i\n',e_conn(i,:)-1);
    end

    fprintf(out_vtk,'\nCELL_TYPES %d\n',n_elements);
    for i=1:n_elements
      fprintf(out_vtk,'22\n');
    end
 
  elseif ( nel_dof==4 )                           % Test for linear quads
    fprintf(out_vtk,'CELLS %7d %7d\n',n_elements,5*n_elements);

    for i=1:n_elements
      fprintf(out_vtk,'4  %6i  %6i  %6i  %6i\n',e_conn(i,:)-1);
    end

    fprintf(out_vtk,'\nCELL_TYPES %d\n',n_elements);
    for i=1:n_elements
      fprintf(out_vtk,'9\n');
    end
    
  elseif ( nel_dof==9 )                           % Test for quadratic quads
    fprintf(out_vtk,'CELLS %7d %7d\n',n_elements,10*n_elements);

    for i=1:n_elements
      fprintf(out_vtk,'9  %6i  %6i  %6i  %6i  %6i  %6i  %6i  %6i  %6i\n',e_conn(i,:)-1);
    end

    fprintf(out_vtk,'\nCELL_TYPES %d\n',n_elements);
    for i=1:n_elements
      fprintf(out_vtk,'28\n');
    end
 
  else
    error('twod_to_vtk: element type not currently implemented\n')
  end

  
  % Write out variables
  %-----------------------------------------------------------------------------
  fprintf(out_vtk,'\nPOINT_DATA %d\n',n_nodes);

  for j=1:n_scalars
    fprintf(out_vtk,'\nSCALARS %s float 1\n',labels{j});
    fprintf(out_vtk,'LOOKUP_TABLE default\n');
    for i=1:n_nodes
      fprintf(out_vtk,' %14.6e\n', scalars(i,j));
    end
  end
  
  for j=1:n_vectors/2
    fprintf(out_vtk,'VECTORS %s float\n',labels{n_scalars+j});
    for i=1:n_nodes
      fprintf(out_vtk,' %14.6e %14.6e %3.1e\n', vectors(i,2*(j-1)+1), ...
                                                vectors(i,2*(j-1)+2), ...
                                                0);
    end
  end
  
  fclose(out_vtk);
end
