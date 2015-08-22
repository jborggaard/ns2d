function [x,e_conn,boundary] = read_msh(mesh_file,surface_ids) 
%-------------------------------------------------------------------------------
%
%  READ_MSH:  A MATLAB function that reads in a Gmsh *.msh file and 
%             extracts nodes, connectivity as well as boundary nodes
%             (optionally associated with surface IDs in the *.geo file)
%
%             Currently assumes a mesh with either:
%                        linear, quadratic, or cubic triangular elements
%                        linear or quadratic tetrahedral elements
%
%             Works with MeshFormat 2.2 ASCII files
%
%  Version: 1.4a
%
%  Usage:  [x,e_conn,boundary_nodes] = read_msh(mesh_file,surface_ids)
%
%  Variables:
%                mesh_file           Gmsh filename (typically *.msh)
%                surface_ids         a vector of surface numbers for 
%                                    identifying associated boundary nodes
%                                    (optional)
%                                    if no surface_ids are specified, then
%                                    physical groups are used to provide
%                                    surfaces and (if applicable) element
%                                    connectivities.
%                x
%                e_conn
%                boundary(k).nodes   a list of all nodes on the
%                                    surface, surface_ids(k)
%                boundary(k).e_conn  connectivities of boundary elements
%                boundary(k).label   physical surface label
%                boundary(#).nodes   a list of all nodes on the boundary
%                                    # is the last number (use 'end')
%
%  Author: Jeff Borggaard, 2013
%
%-------------------------------------------------------------------------------
%  Need to preallocate storage for elem arrays
%  Need to take advantage of clustering of surface element information
%   (physical groups are written together, current implementation leads to 
%   many unnecessary if statements)

  if (nargin==0)
    error('READ_MSH:  A filename is required');
  end

  fid = fopen(mesh_file);
  
  if (nargin==2)    % we only want boundary information for a subset of surfaces
  
    for i=1:3                                %  skip over the header information
      tline = fgetl(fid);
    end
    tline = fgetl(fid);
    
    if ( strcmp( tline(1:5), '$Phys' ) )
      n_phys = fscanf(fid,'%d',1);
      for i=1:n_phys
        tline = fgetl(fid);
      end
      tline = fgetl(fid);
    end
  
    while ( ~strcmp( tline(1:6),'$Nodes') )
      tline = fgetl(fid);
    end

    dim        = 0;   % estimate later
    n_phys     = 0;
    no_phys    = 1;
    n_surfaces = length(surface_ids);
    
  else                   % return all of the physical boundary segments/surfaces
    
    for i=1:3                                %  skip over the header information
      tline = fgetl(fid);
    end
    tline = fgetl(fid);
    
    if ( strcmp( tline(1:5), '$Phys' ) )
      n_phys = fscanf(fid,'%d',1);
      no_phys   = 0;
      for i=1:n_phys
        [ group(i).dim    ] = fscanf(fid,'%d',1);
        [ group(i).number ] = fscanf(fid,'%d',1);
        [ group(i).name   ] = fscanf(fid,'%s',1);
        group(i).name = group(i).name(2:end-1);   % strip off the ""
      end
      tline = fgetl(fid);

      %  Calculate the dimension of the "volume" elements (either 2 or 3)
      dim = 0;
      for i=1:n_phys
        dim = max( dim, group(i).dim );
      end

      tline = fgetl(fid);
      tline = fgetl(fid);  % should be '$Nodes...'

    else  % no physical groups, just return all boundary nodes
      dim         = 0;   % estimate later
      n_phys      = 0;
      surface_ids = [];
    end
      
    surf_count  = 0;
    surface_ids = [];

    for i=1:n_phys
      if (group(i).dim == dim-1)        % this is a physical surface on boundary
        surf_count = surf_count + 1;
        surface_ids(surf_count) = group(i).number;
        boundary(surf_count).name = group(i).name;
      end
    end
    
    boundary(surf_count+1).name = 'all';
    
%     for i=1:n_phys
%       disp(group(i).dim)
%       disp(group(i).number)
%       disp(group(i).name)
%     end
 
    n_surfaces = length(surface_ids);  % these are physical surfaces
      
  end
  
  %%----------------------------------------------------------------------------
  %  Read in nodal coordinates
  %-----------------------------------------------------------------------------
  n_nodes = fscanf(fid,'%d',1);
  
  % Read in the first node to estimate the dimension of the problem if it
  % hasn't been identified.
  [temp] = fscanf(fid,'%d %g %g %g',4);
  if (dim==0)
    if (temp(4)==0)
      dim = 2;
    else
      dim = 3;
    end
    fprintf('Reading in a %d dimensional mesh\n',dim)
  end
  
  x = zeros(n_nodes,dim);

  x(1,1:dim) = temp(2:dim+1);
  
  for n=2:n_nodes
    [temp] = fscanf(fid,'%d %g %g %g',4);
    x(n,1:dim) = temp(2:dim+1);
  end

  tline = fgetl(fid);
  tline = fgetl(fid);
  while ( ~strcmp( tline(1:5),'$Elem') )
    tline = fgetl(fid);
  end

  t_elem = fscanf(fid,'%d',1); % maximum number of 0d, 1d, 2d and/or 3d elements
  
  %%----------------------------------------------------------------------------
  %  Create storage for element connectivity and boundary nodes
  %-----------------------------------------------------------------------------
  e_conn = zeros(t_elem,dim+1);
  nelem  = zeros(1,n_surfaces);

  last = n_surfaces+1;
  boundary(last).nodes = zeros(n_nodes,1); % set all boundary nodes to .false.
  for i=1:max(1,n_surfaces)
    out(i).nodes = zeros(n_nodes,1);
  end
   
  tline = fgetl(fid);
  n_elem = 0;
  for n=1:t_elem
    [elm_type] = fscanf(fid,'%d %d',2);
    
    if     (elm_type(2)==1)                                       %  2-node line
      [temp] = fscanf(fid,'%d %d %d %d %d',5);
    
      if ( dim==2 )
        boundary(last).nodes(temp(4:5))=1;
        
        for k=1:n_surfaces
          if ( temp(2+no_phys)==surface_ids(k) )
            out(k).nodes(temp(4:5))=1;
            nelem(k) = nelem(k) + 1;
            boundary(k).elem(nelem(k),:) = temp(4:5);
          end
        end
      end
      
    elseif (elm_type(2)==2)                                   %  3-node triangle
      [temp] = fscanf(fid,'%d %d %d %d %d %d',6);
    
      if ( dim==2 )
        n_elem = n_elem + 1;
        e_conn(n_elem,1:3) = temp(4:6);

      else
        boundary(last).nodes(temp(4:6))=1;
    
        for k=1:n_surfaces
          if (temp(2+no_phys)==surface_ids(k))
            out(k).nodes(temp(4:6))=1;
            nelem(k) = nelem(k) + 1;
            boundary(k).elem(nelem(k),:) = temp(4:6);
          end
        end
      end      

    elseif (elm_type(2)==4)                         %  4-node linear tetrahedron
      [temp] = fscanf(fid,'%d %d %d %d %d %d %d',7);
      n_elem = n_elem + 1;
      e_conn(n_elem,1:4) = temp(4:7);

      
    elseif (elm_type(2)==8)                          %  3-node second order line
      [temp] = fscanf(fid,'%d %d %d %d %d %d',6);
      
      if ( dim==2 )
        boundary(last).nodes(temp(4:6))=1;
    
        for k=1:n_surfaces
          if (temp(2+no_phys)==surface_ids(k)) 
            out(k).nodes(temp(4:6))=1;
            nelem(k) = nelem(k) + 1;
            boundary(k).elem(nelem(k),:) = temp(4:6);
          end
        end
      end
      
    elseif (elm_type(2)==9)                      %  6-node second order triangle
      [temp] = fscanf(fid,'%d %d %d %d %d %d %d %d %d',9);
    
      if ( dim==2 )
        n_elem = n_elem + 1;
        e_conn(n_elem,1:6) = temp(4:9);

      else
        boundary(last).nodes(temp(4:9))=1;
    
        for k=1:n_surfaces
          if (temp(2+no_phys)==surface_ids(k))
            out(k).nodes(temp(4:9))=1;
            nelem(k) = nelem(k) + 1;
            boundary(k).elem(nelem(k),:) = temp(4:9);
          end
        end
      end
    
    elseif (elm_type(2)==11)                  % 10-node second order tetrahedron
      [temp] = fscanf(fid,'%d %d %d %d %d %d %d %d %d %d %d %d %d',13);
      n_elem = n_elem + 1;
      e_conn(n_elem,1:10) = temp(4:13);

    
    elseif (elm_type(2)==15)                                     %  1-node point
      [temp] = fscanf(fid,'%d %d %d %d',4);
      

    elseif (elm_type(2)==21)                      % 10-node third order triangle
      [temp] = fscanf(fid,'%d %d %d %d %d %d %d %d %d %d %d %d %d',13);

      if ( dim==2 )
        n_elem = n_elem + 1;
        e_conn(n_elem,1:10) = temp(4:13);
      else
        warning('cubic tetrahedra haven''t been implemented')
      end


    else
      error('READ_MSH:  Unexpected elm_type = %d\n',elm_type(2))
      
    end
  end

  e_conn = e_conn(1:n_elem,:);
  
  %  Gmsh has an unconventional node ordering for quadratics, fix it...
  if ( size(e_conn,2)==10 && elm_type(2)==11 )
    e_conn(:,8:10) = e_conn(:,[10 9 8]);
  end
  
  index = (1:n_nodes)';
  boundary(last).nodes = index(boundary(last).nodes>0);
  for k=1:n_surfaces
    out(k).nodes = index(out(k).nodes>0);
    boundary(k).nodes = out(k).nodes;
  end
end
