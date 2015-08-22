function [x,e_conn,boundary] = twod_read_msh(mesh_file,surface_ids) 
%%------------------------------------------------------------------------------
%  TWOD_READ_MSH:  A MATLAB function that reads in a Gmsh *.msh file and 
%                  extracts nodes, connectivity as well as boundary nodes
%                  (optionally associated with surface IDs in the *.geo file)
%
%                  Currently assumes a 2D TRIANGULAR mesh.
%                  Does NOT use physical groups
%
%  Usage:  [x,e_conn,boundary_nodes] = twod_read_msh(mesh_file,surface_ids)
%
%  Variables:
%                x
%                e_conn
%                boundary(1).nodes   a list of all nodes on the boundary
%                boundary(k).nodes   k>1, a list of all nodes on the
%                                    surface, surface_ids(k-1)
%
%  Author: Jeff Borggaard, 2011
%
%-------------------------------------------------------------------------------

  if (nargin==0)
    mesh_file = 'conference_room.msh';
    surface_ids = [3071];
  end

  n_surfaces = length(surface_ids);
  
  fid = fopen(mesh_file);

  %%----------------------------------------------------------------------------
  %  Skip over the header information
  %-----------------------------------------------------------------------------
  tline = fgetl(fid);
  while ( ~strcmp( tline(1:6),'$Nodes') )
    tline = fgetl(fid);
  end

  %%----------------------------------------------------------------------------
  %  Read in nodal coordinates
  %-----------------------------------------------------------------------------
  n_nodes = fscanf(fid,'%d',1);
  x = zeros(n_nodes,2);
  for n=1:n_nodes
    [temp] = fscanf(fid,'%d %g %g %g',4);
    x(n,1:2) = temp(2:3);
  end

  tline = fgetl(fid);
  tline = fgetl(fid);
  while ( ~strcmp( tline(1:5),'$Elem') )
    tline = fgetl(fid);
  end

  t_elem = fscanf(fid,'%d',1);  % maximum number of 0d, 1d and 2d elements
  
  %%----------------------------------------------------------------------------
  %  Create storage for element connectivity and boundary nodes
  %-----------------------------------------------------------------------------
  e_conn = zeros(t_elem,3);

  boundary(1).nodes = zeros(n_nodes,1);
  for i=1:max(1,n_surfaces)
    out(i).nodes = zeros(n_nodes,1);
  end
   
  tline = fgetl(fid);
  n_elem = 0;
  for n=1:t_elem
    [elm_type] = fscanf(fid,'%d %d',2);
    if     (elm_type(2)==1)    %  2-node line
      [temp] = fscanf(fid,'%d %d %d %d %d',5);
      %  temp(3) is the surface id, temp(4:5) are the nodes creating the line.
      boundary(1).nodes(temp(4:5))=1;
      
      for k=1:n_surfaces
        if (temp(3)==surface_ids(k))  % get surface nodes if desired
          out(k).nodes(temp(4:5))=1;
        end
      end
    
      
    elseif (elm_type(2)==2)    %  3-node triangle
      [temp] = fscanf(fid,'%d %d %d %d %d %d',6);
      n_elem = n_elem + 1;
      e_conn(n_elem,1:3) = temp(4:6);
    

    elseif (elm_type(2)==4)    %  4-node linear tetrahedron
      [temp] = fscanf(fid,'%d %d %d %d %d %d %d',7);
      % n_elem = n_elem + 1;
      % e_conn(n_elem,1:4) = temp(4:7);

      
    elseif (elm_type(2)==8)    %  3-node second order line
      [temp] = fscanf(fid,'%d %d %d %d %d %d',6);
      %  temp(3) is the surface id, temp(4:6) are the nodes creating the line.
      boundary(1).nodes(temp(4:6))=1;
    
      for k=1:n_surfaces
        if (temp(3)==surface_ids(k))  % get surface nodes if desired
          out(k).nodes(temp(4:6))=1;
        end
      end
    
      
    elseif (elm_type(2)==9)    %  6-node second order triangle
      [temp] = fscanf(fid,'%d %d %d %d %d %d %d %d %d',9);
      n_elem = n_elem + 1;
      e_conn(n_elem,1:6) = temp(4:9);
    
  
    
    elseif (elm_type(2)==11)  % 10-node second order tetrahedron
      [temp] = fscanf(fid,'%d %d %d %d %d %d %d %d %d %d %d %d %d',13);
      % n_elem = n_elem + 1;
      % e_conn(n_elem,1:10) = temp(4:13);

    
    elseif (elm_type(2)==15)   %  1-node point
      [temp] = fscanf(fid,'%d %d %d %d',4);
      
    else
      error('READ_MSH:  Unexpected elm_type = %d\n',elm_type(2))
      
    end
  end

  e_conn = e_conn(1:n_elem,:);
  
  %  Gmsh has an unconventional node ordering for quadratics, fix it...
  if ( size(e_conn,2)==10 )
    e_conn(:,8:10) = e_conn(:,[10 9 8]);
  end
  
  index = (1:n_nodes)';
  boundary(1).nodes = index(boundary(1).nodes>0);
  for k=1:n_surfaces
    out(k).nodes = index(out(k).nodes>0);
    boundary(k+1).nodes = out(k).nodes;
  end
end
