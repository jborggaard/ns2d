function [fx,fy,q] = f_function0_2d(x,~,~,~)
%% -----------------------------------------------------------------------------
%  Forcing function for the 2D Navier-Stokes equations
%     by convention, the input arguments are: x, t, mu, parameters...
%     the function is called with several quadrature points at once
%     so vector conventions or explicit treatment for each point must be ensured
%-------------------------------------------------------------------------------

  npgaus = size(x,1);
  
  fx = zeros(npgaus,1);
  fy = zeros(npgaus,1);
  q  = zeros(npgaus,1);

end


