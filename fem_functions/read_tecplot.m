function [x,e_conn,data,variables] = read_tecplot(tecplot_file) 
%-------------------------------------------------------------------------------
%  READ_TECPLOT:  A MATLAB function that reads in a tecplot file and 
%                 extracts nodes, connectivity as well as other variables.
%
%                 Assumes ONE ZONE and finite element data types:
%                    TRIANGLE, TETRAHEDRON
%                    QUADRILATERAL, HEXAHEDRON
%                    
%  Usage:  [x,e_conn,data,variables] = read_tecplot(tecplot_file)
%
%  Variables:
%                x
%                e_conn
%                data        a n_nodes x n_variables array containing
%                            the finite element data
%                variables   names of the variables
%
%  Version: 1.3
%
%  Author: Jeff Borggaard, 2013
%
%
%-------------------------------------------------------------------------------

  if (nargin==0)
    error('teplot filename must be provided');
  end
  
  fid = fopen(tecplot_file);
  if (fid==-1)
    error('teplot file %s cannot be opened',tecplot_file);
  end

  %  Skip over the title
  %-----------------------------------------------------------------------------
  [~] = fgetl(fid);

  %  Get variable names (parse a comma-separated string)
  %-----------------------------------------------------------------------------
  tline = fgetl(fid);                      % read VARIABLES = ... string
  tline = tline(strfind(tline,'=')+1:end); % remove "VARIABLES ="

  segments = [0 strfind(tline,',') length(tline)+1];  % parse by ","
  n_var = length(segments)-1;
  
  variables = cell(1,n_var);
  variables{1} = 'X';
  variables{2} = 'Y';
  
  threed = 0;
  if ( n_var==2 )  
    % file only contains 2D mesh data, we are done
    
  else
    n=3;
    token = strtok(tline(segments(n)+1:segments(n+1)-1));
    if ( strcmp(token(1),'"') )
      token = token(2:end-1);
    end
    %   here !
    if ( strcmp(token(1),'Z') || strcmp(token(1),'z') )
      threed = 1;
      variables{3} = 'Z';
    else
      variables{3} = token;
    end
    
    if ( n_var==3 && threed)
      % file only contains 3D mesh data, we are done
    else 
      for n=4:n_var
        token = strtok(tline(segments(n)+1:segments(n+1)-1));
        if ( strcmp(token(1),'"') )
          token = token(2:end-1);
        end
        variables{n} = token;
      end
    end
  end
  
  %  Get element information
  %-----------------------------------------------------------------------------
  tline = fgetl(fid);
  [~,tline] = strtok(tline);   % strips out the keyword 'ZONE'
  
  % assume all keywords are separated by commas
  segments = [0 strfind(tline,',') length(tline)+1];
  
  for i=1:length(segments)-1
    substring = tline(segments(i)+1 : segments(i+1)-1);
    equal_position = strfind(substring,'=');
    keyword = strtok(substring(1:equal_position-1));
    value   = strtok(substring(equal_position+1:end));
    
    if ( strcmp(keyword,'N') )
      n_node = eval(value);
    end
    
    if ( strcmp(keyword,'E') )
      n_elem = eval(value);
    end
    
    if ( strcmp(keyword,'ET') || strcmp(keyword,'ZONETYPE') )
      element_type = value;
    end
    
    if ( strcmp(keyword,'F') || strcmp(keyword,'DATAPACKING') )
      file_type = value;
    end
  end
  
  %  Preallocate storage
%   if ( threed )
%     x = zeros(n_node,3);
%   else
%     x = zeros(n_node,2);
%   end
%  
%   if ( n_var-2-threed>0 )
%     data = zeros(n_node,n_var-2-threed);
%   else
%     data = [];
%   end
 
  if ( strcmp(element_type,'FETRIANGLE') )
    e_conn = zeros(n_elem,3);
    n_dof  = 3;
  elseif ( strcmp(element_type,'TETRAHEDRON') )
    e_conn = zeros(n_elem,4);
    n_dof  = 4;
  elseif ( strcmp(element_type,'QUADRILATERAL') )
    e_conn = zeros(n_elem,8);
    n_dof  = 8;
  end
  
  % test for a blank line
  [tline] = fgetl(fid);
  if (isempty(tline))
    header = 4;
  else
    header = 3;
  end
 
  fclose(fid);
 
  %  Re-open file and read in data
  %-----------------------------------------------------------------------------
  fid = fopen(tecplot_file);

  if ( strcmp(file_type,'POINT') )
    columns = n_var;
    s_format = [];
    for i=1:columns
      s_format = [ s_format ' %f64' ];
    end

    C = textscan(fid, s_format, n_node, 'HeaderLines', header);
   
    if ( threed )
      x(:,1) = C{1};  x(:,2) = C{2};  x(:,3) = C{3};
      for n=4:n_var
        data(:,n-3) = C{n};
      end
    else
      x(:,1) = C{1};  x(:,2) = C{2};
      for n=3:n_var
        data(:,n-2) = C{n};
      end
    end

  elseif ( strcmp(file_type,'FEBLOCK') )
    
    n_lines = ceil(n_node/5);
    s_format = ' %f64 %f64 %f64 %f64 %f64 ';
   
    C = textscan(fid, s_format, n_lines, 'HeaderLines', header);
   
    T(:,1) = C{1};
    T(:,2) = C{2};
    T(:,3) = C{3};
    T(:,4) = C{4};
    T(:,5) = C{5};
    x(:,1) = reshape(T',5*n_lines,1);

    C = textscan(fid, s_format, n_lines, 'HeaderLines', 1);
    T(:,1) = C{1};
    T(:,2) = C{2};
    T(:,3) = C{3};
    T(:,4) = C{4};
    T(:,5) = C{5};
    x(:,2) = reshape(T',5*n_lines,1);

    if ( threed )
      C = textscan(fid, s_format, n_lines, 'HeaderLines', 1);
      T(:,1) = C{1};
      T(:,2) = C{2};
      T(:,3) = C{3};
      T(:,4) = C{4};
      T(:,5) = C{5};
      x(:,3) = reshape(T',5*n_lines,1);
    end
    
    x = x(1:n_node,:);
        
    for n=1:n_var-2-threed
      C = textscan(fid, s_format, n_lines, 'HeaderLines', 1);
      T(:,1) = C{1};
      T(:,2) = C{2};
      T(:,3) = C{3};
      T(:,4) = C{4};
      T(:,5) = C{5};
      data(:,n) = reshape(T',5*n_lines,1);
    end
    
    if ( n_var-2-threed )
      data = data(1:n_node,:);
    end
      
   
  end
 
  %   now read connectivity
  %-----------------------------------------------------------------------------
  s_format = [];
  for i=1:n_dof
    s_format = [ s_format ' %d' ];
  end
   
  C = textscan(fid, s_format, n_elem);
   
  for n=1:n_dof
    e_conn(:,n) = C{n};
  end

 
 
   
end
