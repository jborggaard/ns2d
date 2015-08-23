function [fx,fy] = f_function7_2d(X,t,mu,~)
%% -----------------------------------------------------------------------------
%  A two dimensional verification example motivated by the Stommel ocean flow
%  model.  This was used in Hristova, Etienne, Pelletier and Borggaard 2004
%  to verify unsteady Navier-Stokes simulations (and associated sensitivity
%  analysis with respect to R and beta).
%%------------------------------------------------------------------------------

  npgaus = size(X,1);
  fx = zeros(npgaus,1);
  fy = zeros(npgaus,1);

  % define problem parameters (MUST match the forcing function)
  b      =   1;             % nondimensional latitude parameter
  lambda =   1.591549430918954;          % nondimensional longitude parameter
  D      =   3.1831e-5;     % nondimensional depth
  force  =   1.9249e-11;    % nondimensional wind force
  R      =   1.5198e-6;     % nondimensional friction coefficient
  beta   =   0         ;    % nondimensional Coriolis effect parameter

  F      = force*sin(pi*t); % nonsteady forcing generates a time dependent flow
  Fp     = force*pi*cos(pi*t); % time derivative


  % compute derived constants
  a1 = (D*beta)/(2*R);
  a2 = pi/b;
  k1 = -a1+sqrt( a1^2 + a2^2 );
  k2 = -a1-sqrt( a1^2 + a2^2 );
  C1 = (1-exp(k2*lambda))/(exp(k1*lambda)-exp(k2*lambda));
  C2 = (exp(k1*lambda)-1)/(exp(k1*lambda)-exp(k2*lambda));


  for i = 1:npgaus
    x = X(i,1);
    y = X(i,2);

    u   =      (F /R)*cos(a2*y)*( C1     *exp(k1*x) + C2     *exp(k2*x) - 1 );
    ut  =      (Fp/R)*cos(a2*y)*( C1     *exp(k1*x) + C2     *exp(k2*x) - 1 );
    ux  =      (F /R)*cos(a2*y)*( C1*k1  *exp(k1*x) + C2*k2  *exp(k2*x) );
    uxx =      (F /R)*cos(a2*y)*( C1*k1^2*exp(k1*x) + C2*k2^2*exp(k2*x) );

    uy  =-  a2*(F /R)*sin(a2*y)*( C1     *exp(k1*x) + C2     *exp(k2*x) - 1 );
    uyy =-a2^2*(F /R)*cos(a2*y)*( C1     *exp(k1*x) + C2     *exp(k2*x) - 1 );

    v   =-1/a2*(F /R)*sin(a2*y)*( C1*k1  *exp(k1*x) + C2*k2  *exp(k2*x) );
    vt  =-1/a2*(Fp/R)*sin(a2*y)*( C1*k1  *exp(k1*x) + C2*k2  *exp(k2*x) );
    vx  =-1/a2*(F /R)*sin(a2*y)*( C1*k1^2*exp(k1*x) + C2*k2^2*exp(k2*x) );
    vxx =-1/a2*(F /R)*sin(a2*y)*( C1*k1^3*exp(k1*x) + C2*k2^3*exp(k2*x) );

    vy  =-     (F /R)*cos(a2*y)*( C1*k1  *exp(k1*x) + C2*k2  *exp(k2*x) );
    vyy =   a2*(F /R)*sin(a2*y)*( C1*k1  *exp(k1*x) + C2*k2  *exp(k2*x) );

    px  =-2   *(F /R)*cos(a2*y)*( C1*k1^2*exp(k1*x) + C2*k2^2*exp(k2*x) );
    py  = 2*a2*(F /R)*sin(a2*y)*( C1*k1  *exp(k1*x) + C2*k2  *exp(k2*x) );


    fx(i) = ut + u*ux + v*uy + px - mu*(uxx + uyy);
    fy(i) = vt + u*vx + v*vy + py - mu*(vxx + vyy);

  end
end


