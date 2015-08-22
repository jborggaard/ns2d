  function [] = navierstokes_2d_steady(flow_case)
%% -----------------------------------------------------------------------------
%  navierstokes_2d_steady.m - solves the 2D steady-state Navier-Stokes equations
%
%
%                 Features include:
%                     - Taylor-Hood elements
%                     - Penalty Method
%
%                     - addition of efficient sparse matrix implementation
%
%                     - addition of SUPG terms (in weakjac_2d_steady_supg)
%
%  Copyright (c) 2008, Jeff Borggaard, Virginia Tech
%  Version: 1.3
%
%  Usage:    [] = navierstokes_2d_steady(flow_case)
%
%             
%  Variables:  flow_case
%                         = 1,  validation example
%                               defined in problem01_uv.m
%                         = 2,  validation example
%                               defined in problem02_uv.m
%                         = 5,  flow over a Reuleaux shape
%                               defined in problem05_uv.m
%
%  Internal Variables:
%            Flow Parameters:
%              Re          =  Reynolds number  (1/viscosity)
%
%  Solution parameters
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
  addpath( 'exact'            );
  addpath( 'functions'        );
  addpath( 'problems'         );
  
  
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
    epsilon          = 1e-5/mu;
    
    %  Use mex file to perform element integrations
    if_mex           = true;

    %  Define the initial guess from a file  (usually for using continuation)
    ic_file          = 'none';
     
    %  Define Newton solver parameters
    resid_tol        = 1e-4;
    step_tol         = 1e-4;
    max_iterations   = 25;

    backtracking_iterations =  0;     % =0 is pure Newton...
    lambda_1                = 0.9;
    
    %  Define the Geometry
    [x,e_conn,ide_u,ide_p,dir_u] = problem01_uv(101,21);
    f_function = @f_function0_2d;
    
    matlab_file  = 'solutions_mat/ns_2d_steady_01.mat';
%   tecplot_file = 'solutions_plt/ns_2d_steady_01.plt';
    vtk_file     = 'solutions_vtk/ns_2d_steady_01.vtk';

  
  elseif ( flow_case == 2 )
    %%
    Re               = 100;
    mu               = 1/Re;

    %  Define the penalty parameters
    epsilon          = 1e-5/mu;   % penalty parameter
    
    %  Use mex file to perform element integrations
    if_mex           = true;
    
    %  Define the initial guess from a file  (usually for using continuation)
    ic_file          = 'none';
     
    %  Define Newton solver parameters
    resid_tol        = 1e-4;
    step_tol         = 1e-4;
    max_iterations   = 25;

    backtracking_iterations =  0;     % =0 is pure Newton...
    lambda_1                = 0.9;
    
    %  Define the Geometry
    [x,e_conn,ide_u,ide_p,dir_u] = problem02_uv(21,101);
    f_function = @f_function0_2d;

    matlab_file  = 'solutions_mat/ns_2d_steady_02.mat';
%    tecplot_file = 'solutions_plt/ns_2d_steady_02.plt';
    vtk_file     = 'solutions_vtk/ns_2d_steady_02.vtk';

    
  elseif ( flow_case == 5 )
    %%
    Re               = 100;
    mu               = 1/Re;

    %  Define the penalty parameters
    epsilon          = 1e-5/mu;   % penalty parameter
    
    %  Use mex file to perform element integrations
    if_mex           = true;
    
    %  Define the initial guess from a file  (usually for using continuation)
    ic_file          = 'none';
     
    %  Define Newton solver parameters
    resid_tol        = 1e-4;
    step_tol         = 1e-4;
    max_iterations   = 25;

    backtracking_iterations =  0;     % =0 is pure Newton...
    lambda_1                = 0.9;
    
    %  Define the Geometry
    mesh_root = '/Volumes/borggaard/Software/Meshes/reuleauxrcm';
    [x,e_conn,ide_u,ide_p,dir_u] = problem05_uv(mesh_root);
    f_function = @f_function0_2d;
    
    matlab_file  = 'solutions_mat/ns_2d_steady_05.mat';
%    tecplot_file = 'solutions_plt/ns_2d_steady_05.plt';
    vtk_file     = 'solutions_vtk/ns_2d_steady_05.vtk';
    

  else
    error('navierstokes_2d_steady:  flow_case=%d is not defined\n',flow_case)
  
  end
  
  material.mu         = mu;
  material.epsilon    = epsilon;
  material.f_function = f_function;
  material.if_mex     = logical(if_mex);

  [n_nodes   , n_dimensions] = size(x     );
  [n_elements, nel_dof     ] = size(e_conn);
  
  if ( n_dimensions ~= 2 )
    error('navierstokes_2d_steady: expecting a 2d mesh\n')
  end
  
  if ( nel_dof ~= 6 )
    error('navierstokes_2d_steady: assuming a P2-P1 Taylor-Hood element pair\n')
  end

  %% ---------------------------------------------------------------------------
  %  Writeout Problem Parameters
  %%----------------------------------------------------------------------------
  fprintf('navierstokes_2d_steady:\n')
  fprintf('  Running flow_case %d\n',flow_case)
  fprintf('  The mesh has %d elements and %d nodes\n',n_elements,n_nodes)
  ndofu = max(max(ide_u));   ndofp = max(max(ide_p))-ndofu;
  fprintf('  corresponding to %d velocity and %d pressure unknowns\n',...
          ndofu, ndofp)
  fprintf('  The Reynolds number is %g\n',Re)
  if ( if_mex )
    fprintf('  Using mex file for element integrations\n')
  end
  
  %% ---------------------------------------------------------------------------
  %  Initialize Newton Iterations
  %%----------------------------------------------------------------------------
  if ( exist('ic_file','var') && ~strcmp(ic_file,'none') )  % read initial guess
    load(ic_file,'u','v','p')
    
  else              %  set initial guess to zero then impose boundary conditions
    u = zeros(n_nodes,1);
    v = zeros(n_nodes,1);
  
    for n=1:n_nodes       %  adjust initial guess to satisfy boundary conditions
      index = ide_u(n,1);
      if (index<0)
        u(n) = dir_u(-index);
      end

      index = ide_u(n,2);
      if (index<0)
        v(n) = dir_u(-index);
      end
    end
  
    % initialize the pressure
    p = zeros(n_nodes,1);
  end
  
  t0=toc;
  fprintf('time for setup is: %d\n',t0)

  
  %% ---------------------------------------------------------------------------
  %  Begin Newton Iterations
  %-----------------------------------------------------------------------------
  
  tic  
  [J,res] = weakjac_2d_steady( x, e_conn,  ...
                               material,   ...
                               ide_u,ide_p,...
                               u, v, p        );  
  t1=toc;
  fprintf('  time for matrix integration is: %d\n',t1)
  
  norm_resid = norm(res);

  converged = false;
  diverged  = false;
  iteration = 0;

  while (~converged && ~diverged)
  
    iteration = iteration + 1;
    
    fprintf('--------------------------------------------------\n')
    fprintf('! Newton Iteration %d \n',iteration)
    fprintf('--------------------------------------------------\n')
    
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

        
      [J,res] = weakjac_2d_steady( x, e_conn,                    ...
                                   material,                     ...
                                   ide_u,            ide_p,      ...
                                   u_trial, v_trial, p_trial        );
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
%           % stop?     it didn't converge, but accept what we have...
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
    fprintf(' converge1 = %d, and converge2  = %d\n\n',resid_tol, step_tol);
    converged = ( ( norm_resid < resid_tol & norm_step < step_tol ) );
  
    if (~converged)
      diverged  = iteration>max_iterations;
    end
      
  end % Newton iterations

  if ( diverged )
    fprintf('solution diverged...\n');
  end
  
  %% ---------------------------------------------------------------------------
  %  Post Processing Module
  %-----------------------------------------------------------------------------

  for n=1:n_elements
    local_nodes = e_conn(n,:);
    p( local_nodes(4) ) = ( p( local_nodes(1) ) + p( local_nodes(2) ) )/2;
    p( local_nodes(5) ) = ( p( local_nodes(2) ) + p( local_nodes(3) ) )/2;
    p( local_nodes(6) ) = ( p( local_nodes(3) ) + p( local_nodes(1) ) )/2;
  end

  if ( exist('matlab_file','var') )
    fprintf('Writing solution in %s\n',matlab_file)
    save( matlab_file, 'x', 'e_conn', 'ide_u', 'u', 'v', 'p' )
  end
  
  if ( exist('tecplot_file','var') )
    fprintf('Writing solution in %s\n',tecplot_file)
    twod_to_tecplot( tecplot_file, x, e_conn, [u,v,p], 'navier stokes flow' )
  end

  if ( exist('vtk_file','var') )
    fprintf('Writing solution in %s\n',vtk_file)
    twod_to_vtk(vtk_file,x,e_conn,p,[u v],{'pressure','velocity'})
  end

end % function

%% -----------------------------------------------------------------------------
%  Supporting Functions
%-------------------------------------------------------------------------------

% twod_gauss
% twod_shape
% twod_bilinear
% twod_f_int
% twod_plotc
% twod_plotm1
% twod_plotm2
% twod_getfaces

