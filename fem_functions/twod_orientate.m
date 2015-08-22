function [ t ] = twod_orientate( p,t )
%TWOD_ORIENTATE Ensures that a mesh is positively oriented
%   Usage:  [ e_conn ] = twod_orientate( x, e_conn );
%
%   Variables:
%                    x - vertices of tetrahedralization
%               e_conn - on input: initial connectivity
%                      - on output: positively oriented connectivity

%--------------------------------------------------------------------------
% loop over elements and flip edges for positive orientation if necessary

[n_elems,dof] = size(t);

for n_el=1:n_elems
  v1 = p(t(n_el,2),:)-p(t(n_el,1),:);
  v2 = p(t(n_el,3),:)-p(t(n_el,1),:);
  
  area = v1(1)*v2(2) - v1(2)*v2(1);
  
  if (area<eps)
    t(n_el,2:3) = [t(n_el,3) t(n_el,2)];
  end

  if ( abs(area)< 1e-8 )
    error(' problem with triangle %d\n',n_el);
    error('  %d  %d  %d\n',p(t(n_el,1),:));
    error('  %d  %d  %d\n',p(t(n_el,2),:));
    error('  %d  %d  %d\n',p(t(n_el,3),:));
  end

end


end

