  function [] = navierstokes_2d(flow_case)
%% -----------------------------------------------------------------------------
%  navierstokes_2d.m - solves the 2D time-dependent Navier-Stokes equations
%                      using Taylor-Hood finite elements, Crank-Nicolson
%                      integration, and a penalty method formulation.
%
%                      element integrations are optionally implemented using
%                      mex files.
%
%
%  Copyright (c) 2008, Jeff Borggaard, Virginia Tech
%  Version: 1.3
%
%  Usage:    [] = navierstokes_2d(flow_case)
%
%  Variables:  flow_case
%                          = 1  impulsively started flow in a horizontal channel
%                               defined in problem01_uv.m
%                          = 2  impulsively started flow in a vertical channel
%                               defined in problem02_uv.m
%                          = 5  flow past a Reuleaux triangle (a shape with 
%                               the same cross-sectional area at every angle 
%                               to avoid blockage area changes with angle)
%                               defined in problem05_uv.m
%                          = 7  a Stommel flow verification problem described
%                               in Hristova, et al., 2006.
%                               defined in problem07_uv.m
%                          = 8  flow past a circular cylinder, problem described
%                               in Hay, Borggaard, and Pelletier, 2009.
%                               defined in problem08_uv.m
%
%  Internal Variables:
%            Flow Parameters:
%              Re          =  Reynolds number  (1/viscosity)
%
%            Time Integration Parameters:
%              time_vector =  predefined time integration points, e.g.
%                             linspace( t_initial, t_final, n_steps )
%              n_save      =  number of time steps between saves
%              t_step      =  time step (should be consistent with spacing in
%                             time_vector)
%
%              theta       =  Crank-Nicolson parameter (theta=0.5 is optimal)
%
%            Solution parameters
%              epsilon     =  penalty parameter
%
%              if_mex      =  flag to use a mex function to perform element
%                             integrations
%
%              resid_tol   =  Newton convergence/divergence parameters
%              step_tol
%              max_iterations
%
%              backtracking_iterations = Newton globalization parameters
%              lambda_1
%% -----------------------------------------------------------------------------

  %  define path to the fem_functions library
  addpath( 'fem_functions' );
  addpath( 'exact'         );
  addpath( 'functions'     );
  addpath( 'problems'      );
  
  tic
  
  %% ---------------------------------------------------------------------------
  %  Define "Input" Parameters
  %%----------------------------------------------------------------------------
  if (nargin<1)
    flow_case      = 1;
  end
  
  
  if ( flow_case == 1 )
    %% -------------------------------------------------------------------------
    %  FLOW_CASE == 1:   Horizontal channel flow, low Reynolds number
    %                    used for verification
    %%--------------------------------------------------------------------------
    %  Define the flow parameters
    Re               = 100;
    mu               = 1/Re;

    %  Define the penalty parameters
    epsilon          = 1e-4/mu;
    
    %  Use mex file to perform element integrations
    if_mex           = true;
    
    %  Define time integration parameters
    t_initial        = 0;   
    t_final          = 20;
    n_steps          = 5001;
    n_save           = 50;
    time_vector      = linspace(t_initial,t_final,n_steps);
    t_step           = time_vector(2)-time_vector(1);

    theta            = 0.5;
    
    %  Define the initial condition file
    ic_file          = 'none';
    time_varying_bc  = false;
    
    %  Define Newton solver parameters
    resid_tol        = 1e-4;
    step_tol         = 1e-4;
    max_iterations   = 20;

    backtracking_iterations = 2;
    lambda_1                = 0.9;

    %  Define the Geometry
    [x,e_conn,ide_u,ide_p,dir_u] = problem01_uv(51,31);
    f_function = @f_function0_2d;

    %  Define output options
    matlab_file  = 'solutions_mat/ns_2d_01.mat';
    tecplot_file = 'solutions_plt/ns_2d_01.plt';
    vtk_file     = 'solutions_vtk/ns_2d_01.vtk';
    snap_root    = 'solutions_vtk/ns_2d_01_';
    

  elseif ( flow_case == 2 )
    %% -------------------------------------------------------------------------
    %  FLOW_CASE == 2:  Vertical channel flow, low Reynolds number
    %                   used for verification
    %%--------------------------------------------------------------------------
    %  Define the flow parameters
    Re               = 100;
    mu               = 1/Re;

    %  Define the penalty parameters
    epsilon          = 1e-4/mu;
    
    %  Use mex file to perform element integrations
    if_mex           = true;
    
    %  Define the initial condition file
    ic_file          = 'none';
    time_varying_bc  = false;
    
    %  Define time integration parameters
    t_initial        = 0;   
    t_final          = 20;
    n_steps          = 5001;
    n_save           = 50;
    time_vector      = linspace(t_initial,t_final,n_steps);
    t_step           = time_vector(2)-time_vector(1);

    theta            = 0.5;
    
    %  Define Newton solver parameters
    resid_tol        = 1e-4;
    step_tol         = 1e-4;
    max_iterations   = 20;

    backtracking_iterations = 2;
    lambda_1                = 0.9;

    %  Define the Geometry
    [x,e_conn,ide_u,ide_p,dir_u] = problem02_uv(21,101);
    f_function = @f_function0_2d;

    %  Define output options
    matlab_file  = 'solutions_mat/ns_2d_02.mat';
    tecplot_file = 'solutions_plt/ns_2d_02.plt';
    vtk_file     = 'solutions_vtk/ns_2d_02.vtk';
    snap_root    = 'solutions_vtk/ns_2d_02_';


  elseif ( flow_case == 5 )
    %% -------------------------------------------------------------------------
    %  FLOW_CASE == 5:  Flow past a Reuleaux triangle (a shape with the same
    %                   cross-sectional area at every angle to avoid blockage
    %                   area changes with angle)
    %%--------------------------------------------------------------------------
    %  Define the flow parameters
    Re               = 100;
    mu               = 1/Re;
    
    %  Define the penalty parameters
    epsilon          = 5e-4/mu;
    
    %  Use mex file to perform element integrations
    if_mex           = true;
    
    %  Define time integration parameters
    t_initial        = 0;   
    t_final          = 200;
    n_steps          = 1001;
    n_save           = 10;
    time_vector      = linspace(t_initial,t_final,n_steps);
    t_step           = time_vector(2)-time_vector(1);

    theta            = 0.5;
    
    %  Define the initial condition file
    ic_file          = 'none';%'solutions_mat/ns_steady_2d_05.mat';
    time_varying_bc  = false;
    
    %  Define Newton solver parameters
    resid_tol        = 1e-4;
    step_tol         = 1e-4;
    max_iterations   = 3;

    backtracking_iterations = 1;
    lambda_1                = 1;

    mesh_root = 'problems/reuleauxrcm';
    [x,e_conn,ide_u,ide_p,dir_u] = problem05_uv(mesh_root);
    f_function = @f_function0_2d;
    
    matlab_file  = 'solutions_mat/ns_2d_05.mat';
    tecplot_file = 'solutions_plt/ns_2d_05.plt';
    vtk_file     = 'solutions_vtk/ns_2d_05.vtk';
    snap_root    = 'solutions_vtk/ns_2d_05_';


  elseif ( flow_case == 7 )
    %% -------------------------------------------------------------------------
    %  FLOW_CASE == 7:  The unsteady Stommel verification example
    %%--------------------------------------------------------------------------
    %  Define the flow parameters
    Re               = .2;    % otherwise, p should be scaled by mu
    mu               = 1/Re;
    
    %  Define the penalty parameters
    epsilon          = 1e-12/mu;
    
    %  Use mex file to perform element integrations
    if_mex           = false;
    
    %  Define time integration parameters
    t_initial        = 0;   
    t_final          = .25;
    n_steps          = 400;
    n_save           = 10;
    time_vector      = linspace(t_initial,t_final,n_steps+1);
    t_step           = time_vector(2)-time_vector(1);

    theta            = 0.5;
    
    %  Define the initial condition file
    ic_file          = 'exact';
    time_varying_bc  = true;
    
    %  Define Newton solver parameters
    resid_tol        = 1e-12;
    step_tol         = 1e-12;
    max_iterations   = 10;

    backtracking_iterations = 1;
    lambda_1                = 1;

    [x,e_conn,ide_u,ide_p,dir_u] = problem07_uv(65);
    d_function = @exact_solution7;
    e_function = @exact_solution7;
    f_function = @f_function7_2d;

    matlab_file  = 'solutions_mat/ns_2d_07.mat';
    %tecplot_file = 'solutions_plt/ns_2d_07.plt';
    vtk_file     = 'solutions_vtk/ns_2d_07.vtk';
    snap_root    = 'solutions_vtk/ns_2d_07_';

  elseif ( flow_case == 8 )
    %% -------------------------------------------------------------------------
    %  FLOW_CASE == 8:  The unsteady flow past a cylinder
    %%--------------------------------------------------------------------------
    %  Define the flow parameters
    Re               = 100;
    mu               = 1/Re;
    
    %  Define the penalty parameters
    epsilon          = 5e-5/mu;
    
    %  Use mex file to perform element integrations
    if_mex           = true;
    
    %  Define time integration parameters
    t_initial        = 0;   
    t_final          = 5;
    n_steps          = 301;
    n_save           = 1;
    time_vector      = linspace(t_initial,t_final,n_steps);
    t_step           = time_vector(2)-time_vector(1);

    theta            = 0.5;
    
    %  Define the initial condition file
    ic_file          = 'none';
    time_varying_bc  = false;
    
    %  Define Newton solver parameters
    resid_tol        = 1e-5;
    step_tol         = 1e-5;
    max_iterations   = 5;

    backtracking_iterations = 1;
    lambda_1                = 1;

    [x,e_conn,ide_u,ide_p,dir_u] = problem08_uv();
    f_function = @f_function0_2d;

    matlab_file  = 'solutions_mat/ns_2d_08.mat';
    tecplot_file = 'solutions_plt/ns_2d_08.plt';
    vtk_file     = 'solutions_vtk/ns_2d_08.vtk';
    snap_root    = 'solutions_vtk/ns_2d_08_';

  else
    error('navierstokes_2d:  flow_case=%d is not defined\n',flow_case)
  
  end
  
  material.mu         = mu;
  material.epsilon    = epsilon;
  material.f_function = f_function;
  material.if_mex     = logical(if_mex);  
  material.dt         = t_step;
  material.theta      = theta;
  
  [n_nodes   , n_dimensions] = size(x     );
  [n_elements, nel_dof     ] = size(e_conn);
  
  if ( n_dimensions ~= 2 )
    error('navierstokes_2d:  expecting 2d mesh\n')
  end
  
  if ( nel_dof ~= 6 )
    error('navierstokes_2d:  assuming a P2-P1 Taylor-Hood element pair\n')
  end

  if ( ~exist('n_store','var') )
    n_store = 0;
  end

  %% ---------------------------------------------------------------------------
  %  Set Initial Conditions for the Solution
  %%----------------------------------------------------------------------------
  if ( strcmp(ic_file,'none') )         % initialize from zero, set b.c.'s below
    uc = zeros(n_nodes,1); 
    vc = zeros(n_nodes,1); 
    pc = zeros(n_nodes,1);
  
  elseif ( strcmp(ic_file,'exact') ...
           && exist('e_function','var') )    % initialize with analytic function
    [uc,vc,pc] = e_function(x,0);
    
  else            % load a (e.g. steady-state) solution as the initial condition
    load(ic_file,'u','v','p'); 
    uc = u;    %#ok
    vc = v;    %#ok
    pc = p;    %#ok
 
  end
  
  %  Set boundary conditions at time = t_initial
  %  When specified using a (time-varying) function
  if ( exist('d_function','var') )
    for n=1:n_nodes
      if ( ide_u(n,1)<0 && ide_u(n,2)<0 )
        [uc(n),vc(n)] = d_function(x(n,:),time_vector(1));
      elseif ( ide_u(n,1)<0 )
        [uc(n),  ~  ] = d_function(x(n,:),time_vector(1));
      elseif ( ide_u(n,2)<0 )
        [  ~  ,vc(n)] = d_function(x(n,:),time_vector(1));
      end
    end
  else
    for n=1:n_nodes
      index = ide_u(n,:);
    
      if (index(1)<0)
        uc(n) = dir_u(-index(1));
      end

      if (index(2)<0)
        vc(n) = dir_u(-index(2));
      end
    end
  end
  
%  snapshot_file = strcat( strcat( snap_root,int2str(n_store),'.plt'));

%  twod_to_tecplot( snapshot_file, x, e_conn, [uc,vc,pc], int2str(1) );

  t0=toc;
  fprintf('time for setup is: %d\n',t0)

  str_time = sprintf('%04d',n_store);
  snapshot_file_vtk = strcat( strcat( snap_root,str_time,'.vtk'));
  twod_to_vtk    ( snapshot_file_vtk, x, e_conn, pc,[uc,vc], ...
                              {'pressure','velocity'} );
  fprintf('  Successfully wrote solution file: %s\n', snapshot_file_vtk);

  snapshot_file_mat = strcat( strcat( snap_root,str_time,'.mat'));
  save( snapshot_file_mat, 'u','v','p', 'Re', 'ide_u','ide_p','x','e_conn' )

  %% ---------------------------------------------------------------------------
  %  Solver Module:  Begin time integration using the theta method
  %%----------------------------------------------------------------------------
  for nt=2:length(time_vector)
    time = time_vector(nt);
  
    fprintf('--------------------------------------------------\n')
    fprintf('     time: %g -> %g\n',time_vector(nt-1),time)

    %% -------------------------------------------------------------------------
    %  Set initial guess w/ boundary conditions for the next time step
    %%--------------------------------------------------------------------------
    u = uc;   v = vc;   p = pc;  % set initial iterate from previous step

    %  The initial guess, for nt, must satisfy boundary conditions at
    %  time(nt).
    
    if ( time_varying_bc )

      for n=1:n_nodes
        index_u = ide_u(n,1);
        index_v = ide_u(n,2);
        if ( index_u<0 || index_v<0 )
          [u_bc,v_bc] = d_function(x(n,:),time);
          
          if ( index_u<0 )        
            u(n) = u_bc;
          end
          
          if ( index_v<0 )
            v(n) = v_bc;
          end

        end
      
      end
      
    else
      % the fixed boundary conditions at the next time step inherited from 
      % uc and vc

    end      

    %% -------------------------------------------------------------------------
    %  Begin Newton Iterations
    %%--------------------------------------------------------------------------

    %  Set up the first Newton iteration
    
%       [u,v,p] = e_function(x,time_vector(nt));
%       [uc,vc,pc] = e_function(x,time_vector(nt-1));
    [J,res] = weakjac_2d( x, e_conn,        ...
                          material,         ...
                          ide_u,  ide_p,    ...
                          u , v , p ,       ...
                          uc, vc, pc,       ...
                          time_vector(nt-1),...
                          time_vector(nt  )    );

    norm_resid = norm(res);
    
    converged = false;
    diverged  = false;
    iteration = 0;
  
    while (~converged && ~diverged)
      iteration = iteration + 1;
    
      %% -----------------------------------------------------------------------
      %  Compute Newton Step
      %%------------------------------------------------------------------------
      tic;
      newton_step = -J\res;
      t1 = toc;
      fprintf('  time for linear solve is: %d\n',t1)
      
      %% -----------------------------------------------------------------------
      %  Find an acceptable iterate by backtracking
      %%------------------------------------------------------------------------
      valid       = false;
      b_iteration = 0;
      lambda      = lambda_1;

      while (~valid)
        b_iteration = b_iteration + 1;
        u_trial     = u;
        v_trial     = v;
        p_trial     = p;
        
        for n=1:n_nodes
          index = ide_u(n,1);
          if (index>0)
            u_trial(n) = u_trial(n) + lambda*newton_step(index);
          end

          index = ide_u(n,2);
          if (index>0)
            v_trial(n) = v_trial(n) + lambda*newton_step(index);
          end
        end

        for n=1:n_nodes
          index = ide_p(n,1);
          if (index>0)
            p_trial(n) = p_trial(n) + lambda*newton_step(index);
          end
        end

        norm_step   = lambda*norm(newton_step,2);

        [J,res] = weakjac_2d( x, e_conn,                    ...
                              material,                     ...
                              ide_u,            ide_p,      ...
                              u_trial, v_trial, p_trial,    ...
                              uc     , vc     , pc     ,    ...    
                              time_vector(nt-1),            ...
                              time_vector(nt  )                );
        t1=toc-t1;
        fprintf('  time for matrix integration is: %d\n',t1)


        norm_trial = norm(res);
 
        if ( norm_trial < norm_resid*(1-lambda*10^-8) ) 
          % accept the step
          fprintf('  Step is accepted: %g < %g\n',norm_trial,norm_resid)
          valid      = true;
          u          = u_trial;
          v          = v_trial;
          p          = p_trial;
          norm_resid = norm_trial;

        else
          % Armijo search (backtracking)
          lambda = lambda / 2;
          if ( b_iteration == backtracking_iterations+1 )
            valid = true; %diverged = 1;
            fprintf(['Armijo step too small, accept current and continue\n'...
                     'residual norm was : ' num2str(norm_resid) '\n'])
%             % stop?     it didn't converge, but accept what we have...
            u = u_trial;
            v = v_trial;
            p = p_trial;
            norm_resid = norm_trial;
          end
        end
          
      end  % found a valid step in the line search

      
      t2=toc;
      fprintf('  time for Newton step is: %d\n\n',t2)

      fprintf(' || res || = %d, and || step || = %d\n',norm_resid, norm_step);
      fprintf(' converge1 = %d, and converge2  = %d\n\n',resid_tol*t_step, step_tol);
      converged = ( ( norm_resid < resid_tol*t_step & norm_step < step_tol ) );
  
      if (~converged)
        diverged  = iteration>max_iterations;
      end
      
    end % Newton iterations

    if ( diverged )
      fprintf('solution diverged...\n');
    end

    uc = u;  vc = v;  pc = p;

    if ( exist('e_function','var') )
      [ue,ve,pe] = e_function(x,time_vector(nt));
      
      if ( ~exist('M','var') )
        M = compute_mass_matrix(x,e_conn);
      end
      
      erru = sqrt( (u-ue)'*M*(u-ue) );
      errv = sqrt( (v-ve)'*M*(v-ve) );
      nrmu = sqrt(     ue'*M*ue     );
      nrmv = sqrt(     ve'*M*ve     );
      fprintf('At t=%g, the u solution norm is %g, ',time_vector(nt),nrmu)
      fprintf('and nodal error is %g\n',erru)
      fprintf('At t=%g, the v solution norm is %g, ',time_vector(nt),nrmv)
      fprintf('and nodal error is %g\n',errv)
      
      for n=1:n_elements
        local_nodes = e_conn(n,:);
        p( local_nodes(4) ) = ( p( local_nodes(1) ) + p( local_nodes(2) ) )/2;
        p( local_nodes(5) ) = ( p( local_nodes(2) ) + p( local_nodes(3) ) )/2;
        p( local_nodes(6) ) = ( p( local_nodes(3) ) + p( local_nodes(1) ) )/2;
      end
      ep = pe-p;
      eu = ue-u;
      ev = ve-v;
      twod_to_vtk('error.vtk',x,e_conn,[p pe ep],[u v ue ve eu ev],...
                  {'pressure_n','pressure_e','error_p',...
                   'velocity_n','velocity_e','error_v'} );
    end
    
    if (mod(nt-1,n_save)==0)
      n_store = n_store + 1; 
      for n=1:n_elements
        local_nodes = e_conn(n,:);
        p( local_nodes(4) ) = ( p( local_nodes(1) ) + p( local_nodes(2) ) )/2;
        p( local_nodes(5) ) = ( p( local_nodes(2) ) + p( local_nodes(3) ) )/2;
        p( local_nodes(6) ) = ( p( local_nodes(3) ) + p( local_nodes(1) ) )/2;
      end

      str_time = sprintf('%04d',n_store);

      if ( exist('tecplot_file','var') )
        snapshot_file_tec = strcat( strcat( snap_root,str_time,'.plt'));
        twod_to_tecplot( snapshot_file_tec, x, e_conn, [u,v,p], int2str(nt) );
        fprintf('  Successfully wrote solution file: %s\n', snapshot_file_tec);
      end
      
      snapshot_file_vtk = strcat( strcat( snap_root,str_time,'.vtk'));
      twod_to_vtk    ( snapshot_file_vtk, x, e_conn, p,[u,v], ...
                                  {'pressure','velocity'} );
      fprintf('  Successfully wrote solution file: %s\n', snapshot_file_vtk);

      snapshot_file_mat = strcat( strcat( snap_root,str_time,'.mat'));
      save( snapshot_file_mat, 'u','v','p' )
      fprintf('  Successfully wrote solution file: %s\n', snapshot_file_mat);
      
      save fem_checkpoint.mat u v p x e_conn nt time_vector
    end
    
  end % time integration loop

  %% ---------------------------------------------------------------------------
  %  Write out the final solution
  %%----------------------------------------------------------------------------

  for n=1:n_elements
    local_nodes = e_conn(n,:);
    p( local_nodes(4) ) = ( p( local_nodes(1) ) + p( local_nodes(2) ) )/2; %#ok
    p( local_nodes(5) ) = ( p( local_nodes(2) ) + p( local_nodes(3) ) )/2; %#ok
    p( local_nodes(6) ) = ( p( local_nodes(3) ) + p( local_nodes(1) ) )/2; %#ok
  end

  if ( exist('matlab_file','var') )
    fprintf('Writing solution in %s\n',matlab_file)
    save(matlab_file, 'x', 'e_conn', 'u', 'v', 'p')
  end

  if ( exist('tecplot_file','var') )
    fprintf('Writing solution in %s\n',tecplot_file)
    twod_to_tecplot( tecplot_file, x, e_conn, [u,v,p], 'navier stokes flow' )
  end

  if ( exist('vtk_file','var') )
    fprintf('Writing solution in %s\n',vtk_file)
    twod_to_vtk(vtk_file,x,e_conn,p,[u v],{'pressure','velocity'})
    if ( exist('e_function','var') )
      twod_to_vtk('error.vtk',x,e_conn,[p, pe],[u, v, ue, ve, (ue-u), (ve-v)])
    end
  end

end % function

%% -----------------------------------------------------------------------------
%  List of Supporting Functions
%%------------------------------------------------------------------------------

% twod_gauss
% twod_shape
% twod_bilinear
% twod_f_int
