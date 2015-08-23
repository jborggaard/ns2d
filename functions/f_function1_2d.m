function [fx,fy,q] = f_function1_2d(x,~,mu,kappa,beta)
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
  
  k     = kappa;
  expan = beta;

  U0   = 1;
  dens = 1;
  c_p  = 1;
  dT   = 1;
  T0   = 0;
  g    =-1;

  alpha = sqrt(U0*dens/mu);

  n = 1/3;
  beta2 = alpha*(mu*c_p/k)^n;

  for i = 1:npgaus                                                          

    eta   = alpha*x(i,2)/sqrt(x(i,1));
    expo1 = exp(-eta);
    fctx1 = alpha/sqrt(x(i,1));

    gamma = beta2*x(i,2)/sqrt(x(i,1));
    expo2 = exp(-gamma);
    fctx2 = beta2/sqrt(x(i,1));


    u = U0*(1 - expo1);
    v = U0*(1-expo1*(1+eta))/(2*alpha*sqrt(x(i,1)));
    T = dT*expo2+T0;

    ux  = -U0*0.5*eta*expo1/x(i,1);
    uy  = U0*fctx1*expo1;
    uxx = U0*expo1*(3*eta-eta^2)/(4*x(i,1)^2);
    uyy = -U0*fctx1^2*expo1;

    vx  = U0*(expo1*(1+eta-eta^2)-1)*x(i,1)^(-1.5)/(4*alpha);
    vy  = -ux;
    
    vxx = U0*(expo1*(-3-3*eta+6*eta^2-eta^3) + 3)...
                *x(i,1)^(-2.5)/(8*alpha);
    vyy = U0*fctx1*expo1*(1-eta)/(2*x(i,1));

    px = 1;
    py = 0;

    fx(i) = dens*(u*ux + v*uy) + px - mu*(uxx+uyy) + dens*g*expan*(T-T0);
    fy(i) = dens*(u*vx + v*vy) + py - mu*(vxx+vyy);

    Tx = 0.5*dT*gamma*expo2/x(i,1);
    Ty = -dT*fctx2*expo2;

    Txx = -dT*expo2*(3*gamma-gamma^2)/(4*x(i,1)^2);
    Tyy = dT*fctx2^2*expo2;

    q(i) = dens*c_p*(u*Tx + v*Ty) - k*(Txx+Tyy);
  end
end


