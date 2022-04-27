%%%
%%% ACC_params
%%%
%%% Structure to contain model parameters
%%%
classdef ACC_params 
%  
  properties
    N
    Lx
    Ly
    H1
    H2
    Hb
    Wb
    Xb
    f
    g
    rg
    rho0
    tau
    Taabw
    rb
    nu1
    nu2
    K
    Ky
    beta 
 end
%  
  methods
%    
    function obj = ACC_params()      
      obj.N=512; % number of fourier modes
      obj.Lx=3200000; % zonal domain size
      obj.Ly=1600000; % meridional domain size
      obj.H1=1500; % mean upper layer thickness
      obj.H2=2500; % mean bottom layer thickness
      obj.Hb=1000; % ridge height
      obj.Wb=150000; % ridge width
      obj.Xb=0; % ridge longitude
      obj.f=-1e-4; % coriolis parameter
      obj.g=9.81; % reduced gravity
      obj.rg=.01; % reduced gravity
      obj.rho0=1000; % density
      obj.tau = 0.1; % wind stress
      obj.Taabw = 0; % AABW export
      obj.rb=2e-4; % drag velocity
      obj.nu1=1000; % Eddy viscosity
      obj.nu2=1000; % Eddy viscosity
      obj.K=400; % Eddy diffusion
      obj.Ky=0; %%% Meridional eddy diffusivity
      obj.beta=1.5e-11; % coriolis parameter gradient
    end
%        
  end
end