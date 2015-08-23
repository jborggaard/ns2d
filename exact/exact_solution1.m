function [u,v,p,T,  ux,uy, vx,vy, px,py, Tx,Ty] = exact_solution1(x)
%-------------------------------------------------------------------------------
%  Matches analytic solution in AIAA 99-3625 and 2000-4732 (w/ constant mu)
%-------------------------------------------------------------------------------
  npgaus = size(x,1);
  
  p = zeros(npgaus,1);  px = p;  py = p;
  u = zeros(npgaus,1);  ux = u;  uy = u;
  v = zeros(npgaus,1);  vx = v;  vy = v;
  T = zeros(npgaus,1);  Tx = T;  Ty = T;
  
  Re = 100;
  Pr = 0.7;

  U0 = 1;
  dens = 1;
  mu = 1/Re;
  c_p = 1;
  k = 1/(Re*Pr);
  dT = 1;
  T0 = 0;

  alpha = sqrt(U0*dens/mu);

  n = 1/3;
  beta2 = alpha*(mu*c_p/k)^n;

  for i = 1:npgaus                                                         
                                                                                
    xg = x(i,1);
    yg = x(i,2);                                                            

    eta = alpha*yg/sqrt(xg);
    expo1 = exp(-eta);
    fctx1 = alpha/sqrt(xg);

    gamma2 = beta2*yg/sqrt(xg);
    expo2 = exp(-gamma2);
    fctx2 = beta2/sqrt(xg);

    p(i)  =  xg-0.1;
    px(i) =  1;
    py(i) =  0;

    u(i)  =  U0*(1 - expo1);
    v(i)  =  U0*(1-expo1*(1+eta))/(2*alpha*sqrt(xg));
    ux(i) = -U0*0.5*eta*expo1/xg;
    uy(i) =  U0*fctx1*expo1;
         
    vx(i) =  U0*(expo1*(1+eta-eta^2)-1)*xg^(-1.5) /(4*alpha);
    vy(i) = -ux(i);
                                          

    T (i) =  dT*expo2+T0;
    Tx(i) =  0.5*dT*gamma2*expo2/xg;
    Ty(i) = -dT*fctx2*expo2;
  end
  
  
end % function
