%%%
%%% ACC_optimizer.m
%%%
%%% Solves the mean flow/standing wave equations for a given collection of
%%% model parameters, 'params'.
%%%
function U =  ACC_optimizer (params) 

  % extract model parameters
  N=params.N; % number of fourier modes
  Lx=params.Lx; % zonal domain size
  Ly=params.Ly; % meridional domain size
  H1=params.H1; % mean upper layer thickness
  H2=params.H2; % mean bottom layer thickness
  Hb=params.Hb; % ridge height
  Wb=params.Wb; % ridge width
  Xb=params.Xb; % ridge longitude
  f=params.f; % coriolois parameter
  tau=params.tau; % wind stress
  rho0=params.rho0; % density
  rb=params.rb; % drag velocity    
  K=params.K; % Eddy diffusion
  beta=params.beta; % coriolis parameter gradient   
  rg = params.rg; % reduced gravity
  f = params.f; % Coriolis parameter
  H1 = params.H1; % upper layer thickness
  H2 = params.H2; % lower layer thickness
  
  %%% 1 over squared deformation radius for calculation
  Ld1=sqrt(rg*H1)/abs(f); %%% upper layer deformation radius
  Ld2=sqrt(rg*H2)/abs(f); %%% lower layer deformation radius
  Ld1sq=1/Ld1^2; 
  Ld2sq=1/Ld2^2;

  % grids in real and spectral space    
  [x,k,etab,etabhat] = gen_grids (N,Lx,Hb,Xb,Wb);
    
  %%% Call lsqnonlin to optimize solution for U
  U0=[0.05 0.01];
  Umin = [-1 -1];
  Umax = [1 1];
  myoptions = optimoptions(@lsqnonlin, ... 
      'StepTolerance',1e-32, ...
      'FunctionTolerance',1e-18, ...
      'OptimalityTolerance',1e-18, ...
      'MaxFunctionEvaluations',10000, ...
      'Algorithm','trust-region-reflective', ... % 'levenberg-marquardt', ...
      'Display','none');
  [U,resnorm] = lsqnonlin(@fun,U0,Umin,Umax,myoptions);

  function res = fun(U)
       
    % generate Fourier coefficients for standing wave equations
    [c1,c2,c3,c4,c5,z1,z2] = calc_sw_coeffs(U,k,params);

    % compute residuals from mean flow equations
    tfs=-f*real(sum(1i*k.*abs(etabhat).^2.*z2,2));
    ifs=H2*Ld2sq*real(sum(1i*k.*abs(z2).^2.*abs(etabhat).^2.*conj(z1),2));
    res=[tau/rho0 - rb*U(2) - tfs, ...
         -rb*U(2)+H2*K*Ld2sq*(U(1)-U(2))+ifs-tfs];
  end

end