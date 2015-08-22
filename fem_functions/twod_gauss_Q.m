function [r,s,w] = twod_gauss_Q(rule)
%-------------------------------------------------------------------------------
%  twod_gauss.m - calculate Gauss integration points for quadrilateral 
%                 elements.  Tensor products of 1D rules are constructed.
%
%  Copyright (c) 2011, Jeff Borggaard, Virginia Tech
%  Version: 1.0
%
%  Usage:    [r,s,w] = twod_gauss_Q(rule)
%
%  Variables:     rule
%                        Number of Gauss points (in the 1d rule):
%                 r
%                        xi coordinate of Gauss points
%                 s
%                        eta coordinate of Gauss points
%                 w
%                        Gauss weights corresponding to (r,s)
%-------------------------------------------------------------------------------

  [r1,w1] = oned_gauss(rule);
  
  r = reshape(ones(rule,1)*r1',rule^2,1);
  s = reshape(r1*ones(1,rule) ,rule^2,1);
  
  w = reshape(w1*w1',rule^2,1);
  
end
