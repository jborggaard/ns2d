function [] = twod_to_pos(filename,x,e_conn,u)
%------------------------------------------------------------------------------
%  twod_to_pos  - function that accepts 2D finite element data and writes
%                 a post-processing file for Gmsh (although this is an
%                 older format, it is used by the adaptive meshing capabilities
%                 in Gmsh).
%
%  Usage:
%         [] = twod_to_pos(filename,x,e_conn,u)
%
%  Variables:     filename
%                      name for the output *.pos file
%                 x
%                      spatial coordinates  (nodes x 3)
%                 e_conn
%                      element connectivity matrix  (elements x (either 3 or 6))
%                 u
%                      field variables
%------------------------------------------------------------------------------
  [n_node,n_var] = size(u);
  [n_elem,n_dof] = size(e_conn);
  
  if (n_dof~=3)
    fprintf('twod_to_pos: Warning, only linear elements are supported\n\n');
  end

  if (n_var>1)
    fprintf('twod_to_pos: Warning, only the first column of the variable is used\n\n')
  end
  
  [fid] = fopen(filename,'w');

  % Header information
  fprintf(fid,'View "background mesh" {\n');
  
  for i=1:n_elem
    nn = e_conn(i,1:3);
    fprintf(fid,'ST(%11.5g, %11.5g, %11.5g,',x(nn(1),1),x(nn(1),2),0);
    fprintf(fid,'   %11.5g, %11.5g, %11.5g,',x(nn(2),1),x(nn(2),2),0);
    fprintf(fid,'   %11.5g, %11.5g, %11.5g)',x(nn(3),1),x(nn(3),2),0);
    fprintf(fid,'{%11.5g,%11.5g,%11.5g};\n',u(nn(1),1),u(nn(2),1),u(nn(3),1));
  end
  fprintf(fid,'};\n');
  
  fclose(fid);
end % function twod_to_pos
