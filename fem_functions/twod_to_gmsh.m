function [] = twod_to_gmsh(filename,x,e_conn,u,title)
%------------------------------------------------------------------------------
%  twod_to_gmsh - function that accepts 2D finite element data and writes
%                 an input file to display solutions in Gmsh (assumes the
%                 geometry is already available as a Gmsh file).
%
%  Usage:
%         [] = twod_to_gmsh(filename,x,e_conn,u,title)
%
%  Variables:     x
%                      spatial coordinates  (nodes x 3)
%                 e_conn
%                      element connectivity matrix  (elements x (either 3 or 6))
%                 u
%                      field variables
%                 filename
%                      name for the output tecplot file
%                 title
%                      title string required by tecplot
%------------------------------------------------------------------------------
  [n_node,n_var] = size(u);
  [n_elem,n_dof] = size(e_conn);
  
  if (n_dof~=3 & n_dof~=6)
    error('twod_to_gmsh: element type not implemented\n')
  end

  [fid] = fopen(filename,'w');

  % Header information
  fprintf(fid,'$MeshFormat\n');
  fprintf(fid,'2.2 0 8\n');
  fprintf(fid,'$EndMeshFormat\n');
  
  % Node Coordinates
  fprintf(fid,'$Nodes\n');
  fprintf(fid,'%d\n',n_node);
  for i=1:n_node
    fprintf(fid,'%d %11.5g %11.5g 0\n',i,x(i,1),x(i,2));
  end
  fprintf(fid,'$EndNodes\n');
  
  % Element Connectivity
  fprintf(fid,'$Elements\n');
  fprintf(fid,'%d\n',n_elem);
  if (n_dof==3)
    for i=1:n_elem
      fprintf(fid,'%d 2 0 %d %d %d\n',i,e_conn(i,:));
    end
  else
    for i=1:n_elem
      fprintf(fid,'%d 9 0 %d %d %d %d %d %d\n',i,e_conn(i,:));
    end
  end
  fprintf(fid,'$EndElements\n');

  % Nodal Data
  fprintf(fid,'$NodeData\n');
  fprintf(fid,'1\n');
  fprintf(fid,strcat(strcat('Title="',title),'"\n'));
  fprintf(fid,'1\n');
  fprintf(fid,'0.0\n');
  fprintf(fid,'3\n');
  fprintf(fid,'0\n');  % timestep number
  if (n_var==1 | n_var==3 | n_var==9) % no zero padding required
    fprintf(fid,'%d\n',n_var);
    n_varp = n_var;
  elseif (n_var==2)
    fprintf(fid,'%d\n',3);
    n_varp = 3;
  elseif (n_var<9)
    fprintf(fid,'%d\n',9);
    n_varp = 9;
  else
    error('twod_to_gmsh: no more than 9 values per node can be processed\n')
  end
    
  fprintf(fid,'%d\n',n_node);

  print_str = '%d ';
  for i=1:n_varp
    print_str = strcat( print_str, ' %11.5g ' );
  end
  print_str = strcat( print_str, '\n' );

  if (n_var==n_varp)
    for i=1:n_node
      fprintf(fid,print_str,i,u(i,:));
    end
  else
    padding = zeros(1,n_varp-n_var);
    for i=1:n_node
      fprintf(fid,print_str,i,[u(i,:),padding]);
    end
  end
  fprintf(fid,'$EndNodeData\n');
  
  fclose(fid);
  % end function threed_to_gmsh
end
