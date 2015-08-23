function [u , v , p , ux , uy , vx , vy , px , py ,...
          su, sv, sp, sux, suy, svx, svy, spx, sp_y]   = exact_solution7(X,t,R,beta)
%% -----------------------------------------------------------------------------
%  A two dimensional verification example motivated by the Stommel ocean flow
%  model.  This was used in Hristova, Etienne, Pelletier and Borggaard 2004
%  to verify unsteady Navier-Stokes simulations (and associated sensitivity
%  analysis with respect to R and beta).
%%------------------------------------------------------------------------------

  % define problem parameters (MUST match the forcing function)
  b      =   1;                 % nondimensional latitude parameter
  lambda =   1.591549430918954; % nondimensional longitude parameter
  D      =   3.1831e-5;         % nondimensional depth
  force  =   1.9249e-11;        % nondimensional wind force

  if ( nargin<3 )
    R      =   1.5198e-6;     % nondimensional friction coefficient
    beta   =   0         ;    % nondimensional Coriolis effect parameter
  end

  F      = force*sin(pi*t); % nonsteady forcing generates a time dependent flow


  % compute derived constants
  a1 = (D*beta)/(2*R);
  a2 = pi/b;
  k1 = -a1+sqrt( a1^2 + a2^2 );
  k2 = -a1-sqrt( a1^2 + a2^2 );
  C1 = (1-exp(k2*lambda))/(exp(k1*lambda)-exp(k2*lambda));
  C2 = (exp(k1*lambda)-1)/(exp(k1*lambda)-exp(k2*lambda));

  
  np = size(X,1);

  % pre-allocate storage
  u  = zeros(np,1);  ux  = u;   uy  = u;
  v  = zeros(np,1);  vx  = v;   vy  = v;
  p  = zeros(np,1);  px  = p;   py  = p;

  su = zeros(np,2);  sux = su;  suy = su;
  sv = zeros(np,2);  svx = sv;  svy = sv;
  sp = zeros(np,2);  spx = sp;  sp_y = sp;


  % evaluate the analytic solution at the given points
  x = X(:,1);
  y = X(:,2);
  
  u  =      (F/R)*cos(a2*y).*( C1     *exp(k1*x) + C2     *exp(k2*x) - 1 );
  ux =      (F/R)*cos(a2*y).*( C1*k1  *exp(k1*x) + C2*k2  *exp(k2*x)     );
  uy =-  a2*(F/R)*sin(a2*y).*( C1     *exp(k1*x) + C2     *exp(k2*x) - 1 );

  v  =-1/a2*(F/R)*sin(a2*y).*( C1*k1  *exp(k1*x) + C2*k2  *exp(k2*x)     );
  vx =-1/a2*(F/R)*sin(a2*y).*( C1*k1^2*exp(k1*x) + C2*k2^2*exp(k2*x)     );
  vy =-     (F/R)*cos(a2*y).*( C1*k1  *exp(k1*x) + C2*k2  *exp(k2*x)     );

  p  =-2   *(F/R)*cos(a2*y).*( C1*k1  *exp(k1*x) + C2*k2  *exp(k2*x)     );
  px =-2   *(F/R)*cos(a2*y).*( C1*k1^2*exp(k1*x) + C2*k2^2*exp(k2*x)     );
  py = 2*a2*(F/R)*sin(a2*y).*( C1*k1  *exp(k1*x) + C2*k2  *exp(k2*x)     );
  
%   for i = 1:np
%     x = X(i,1);
%     y = X(i,2);
% 
%     u  (i) =      (F/R)*cos(a2*y)*( C1     *exp(k1*x) + C2     *exp(k2*x) - 1 );
%     ux (i) =      (F/R)*cos(a2*y)*( C1*k1  *exp(k1*x) + C2*k2  *exp(k2*x)     );
%     uy (i) =-  a2*(F/R)*sin(a2*y)*( C1     *exp(k1*x) + C2     *exp(k2*x) - 1 );
% 
%     v  (i) =-1/a2*(F/R)*sin(a2*y)*( C1*k1  *exp(k1*x) + C2*k2  *exp(k2*x)     );
%     vx (i) =-1/a2*(F/R)*sin(a2*y)*( C1*k1^2*exp(k1*x) + C2*k2^2*exp(k2*x)     );
%     vy (i) =-     (F/R)*cos(a2*y)*( C1*k1  *exp(k1*x) + C2*k2  *exp(k2*x)     );
% 
%     p  (i) =-2   *(F/R)*cos(a2*y)*( C1*k1  *exp(k1*x) + C2*k2  *exp(k2*x)     );
%     px (i) =-2   *(F/R)*cos(a2*y)*( C1*k1^2*exp(k1*x) + C2*k2^2*exp(k2*x)     );
%     py (i) = 2*a2*(F/R)*sin(a2*y)*( C1*k1  *exp(k1*x) + C2*k2  *exp(k2*x)     );
%   end

  if ( nargout>9 )
    % sensitivity with respect to R
    a1p = -a1/R; %(D*beta)/(2*R^2);
    a2p = 0;
    k1p = -a1p + a1*a1p/sqrt( a1^2 + a2^2 );
    k2p = -a1p - a1*a1p/sqrt( a1^2 + a2^2 );
    C1p = ((lambda*k1p*exp(lambda*k1) ...
        - lambda*k2p*exp(lambda*k2))*(exp(lambda*k2) - 1))/   ...
                    (exp(lambda*k1) - exp(lambda*k2))^2 -    ...
         (lambda*k2p*exp(lambda*k2))/(exp(lambda*k1) - exp(lambda*k2));
    C2p = (lambda*k1p*exp(lambda*k1))/(exp(lambda*k1) - exp(lambda*k2)) - ...
         ((lambda*k1p*exp(lambda*k1) - ...
         lambda*k2p*exp(lambda*k2))*(exp(lambda*k1) - 1))/(exp(lambda*k1) - ...
         exp(lambda*k2))^2;
     
    x = X(:,1);
    y = X(:,2);

    su  = (F*cos(y*a2).*(exp(x*k1)*C1p + exp(x*k2)*C2p ...
                    + x.*exp(x*k1)*k1p*C1 + x.*exp(x*k2)*k2p*C2))/R ...
        - (F*cos(y*a2).*(exp(x*k1)*C1 + exp(x*k2)*C2 - 1))/R^2 ...
        - (F*a2p*y.*sin(y*a2).*(exp(x*k1)*C1 + exp(x*k2)*C2 - 1))/R;
    sux = zeros(np,1);
    suy = zeros(np,1);

    sv  = (F*sin(y*a2).*(exp(x*k1)*k1*C1 + exp(x*k2)*k2*C2))/(R^2*a2) ...
        - (F*sin(y*a2).*(exp(x*k1)*k1*C1p + exp(x*k1)*k1p*C1 ...
                      + exp(x*k2)*k2*C2p + exp(x*k2)*k2p*C2 ...
                      + x.*exp(x*k1)*k1*k1p*C1 + x.*exp(x*k2)*k2*k2p*C2))/(R*a2)...
          + (F*sin(y*a2).*(exp(x*k1)*k1*C1 + exp(x*k2)*k2*C2)*a2p)/(R*a2^2) ...
          - (F*y.*cos(y*a2).*(exp(x*k1)*k1*C1 + exp(x*k2)*k2*C2)*a2p)/(R*a2);
    svx = zeros(np,1);
    svy = zeros(np,1);

    sp  = (2*F*cos(y*a2).*(exp(x*k1)*k1*C1 + exp(x*k2)*k2*C2))/R^2 ...
               - (2*F*cos(y*a2).*(exp(x*k1)*k1*C1p + exp(x*k1)*k1p*C1 ...
               + exp(x*k2)*k2*C2p + exp(x*k2)*k2p*C2 ...
               + x.*exp(x*k1)*k1*k1p*C1 + x.*exp(x*k2)*k2*k2p*C2))/R ...
               + (2*F*y.*sin(y*a2).*(exp(x*k1)*k1*C1 + exp(x*k2)*k2*C2)*a2p)/R;
    spx = zeros(np,1);
    sp_y = zeros(np,1);
%     for i=1:np
%       x = X(i,1);
%       y = X(i,2);
% 
%       su (i,1) = (F*cos(y*a2)*(exp(x*k1)*C1p + exp(x*k2)*C2p ...
%                       + x*exp(x*k1)*k1p*C1 + x*exp(x*k2)*k2p*C2))/R ...
%           - (F*cos(y*a2)*(exp(x*k1)*C1 + exp(x*k2)*C2 - 1))/R^2 ...
%           - (F*y*sin(y*a2)*a2p*(exp(x*k1)*C1 + exp(x*k2)*C2 - 1))/R;
%       sux(i,1) = 0;
%       suy(i,1) = 0;
% 
%       sv (i,1) = (F*sin(y*a2)*(exp(x*k1)*k1*C1 + exp(x*k2)*k2*C2))/(R^2*a2) ...
%           - (F*sin(y*a2)*(exp(x*k1)*k1*C1p + exp(x*k1)*k1p*C1 ...
%                         + exp(x*k2)*k2*C2p + exp(x*k2)*k2p*C2 ...
%                         + x*exp(x*k1)*k1*k1p*C1 + x*exp(x*k2)*k2*k2p*C2))/(R*a2)...
%           + (F*sin(y*a2)*(exp(x*k1)*k1*C1 + exp(x*k2)*k2*C2)*a2p)/(R*a2^2) ...
%           - (F*y*cos(y*a2)*(exp(x*k1)*k1*C1 + exp(x*k2)*k2*C2)*a2p)/(R*a2);
%       svx(i,1) = 0;
%       svy(i,1) = 0;
% 
%       sp (i,1) = (2*F*cos(y*a2)*(exp(x*k1)*k1*C1 + exp(x*k2)*k2*C2))/R^2 ...
%                - (2*F*cos(y*a2)*(exp(x*k1)*k1*C1p + exp(x*k1)*k1p*C1 ...
%                + exp(x*k2)*k2*C2p + exp(x*k2)*k2p*C2 ...
%                + x*exp(x*k1)*k1*k1p*C1 + x*exp(x*k2)*k2*k2p*C2))/R ...
%                + (2*F*y*sin(y*a2)*(exp(x*k1)*k1*C1 ...
%                + exp(x*k2)*k2*C2)*a2p)/R;
%       spx(i,1) = 0;
%       sp_y(i,1) = 0;
%     end
    
%     % sensitivity with respect to beta
%       da1_dp   = D/(2*R);
% 
      su (:,2) = zeros(np,1);
      sux(:,2) = zeros(np,1);
      suy(:,2) = zeros(np,1);

      sv (:,2) = zeros(np,1);
      svx(:,2) = zeros(np,1);
      svy(:,2) = zeros(np,1);

      sp (:,2) = zeros(np,1);
      spx(:,2) = zeros(np,1);
      sp_y(:,2) = zeros(np,1);
%     end
  end


end % function
