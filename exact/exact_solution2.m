function [u,v,p,T,  ux,uy, vx,vy, px,py, Tx,Ty] = exact_solution2(X)
%-------------------------------------------------------------------------------
%  Verification example:  borggaard2000soc, AIAA 2000-0563
%% -----------------------------------------------------------------------------
  npgaus = size(X,1);
  
  p = zeros(npgaus,1);  px = p;  py = p;
  u = zeros(npgaus,1);  ux = u;  uy = u;
  v = zeros(npgaus,1);  vx = v;  vy = v;
  T = zeros(npgaus,1);  Tx = T;  Ty = T;
  
  a     = 5000;
  b     =  400;
                                                                    
  gamma = 2;
  beta  = 4;
  d     = 0.01;

%   dens  = 1;
%   mu    = 1/Re;
%   c_p   = 1;
%   k     = 1/(Re*Pr);
%   dT    = 1;
%   T0    = 0;

  for i = 1:npgaus                                                         
    x = X(i,1);  
    y = X(i,2);                                                             

    alpha0 = a*x*x*x*x*y*y;
    alpha1 = gamma*x/( d*(1+beta*y*y) );

    p (i) =  exp(-b*(x*x+y*y));
    px(i) = -2*b*x*p(i);
    py(i) = -2*b*y*p(i);

    T (i) =  exp( alpha1 );
    Tx(i) =  gamma           /(d*(1+beta*y*y)  )*T(i);
    Ty(i) = -2*gamma*x*beta*y/(d*(1+beta*y*y)^2)*T(i);

    u (i) =  4*x*x*    exp( -alpha0 )*( 1 - 2*alpha0              );
    ux(i) =  8*x*      exp( -alpha0 )*( 1 - 8*alpha0 + 4*alpha0^2 );
    uy(i) =  8*a*x^6*y*exp( -alpha0 )*(-3 + 2*alpha0              );

    v (i) = -8*x*y*exp( -alpha0 )*( 1 -  2*alpha0              );
    vx(i) = -8*  y*exp( -alpha0 )*( 1 - 14*alpha0 + 8*alpha0^2 );
    vy(i) = -ux(i);
  end
  
  
end % function
