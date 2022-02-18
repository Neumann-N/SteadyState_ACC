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
    nu
    K
    beta
    Ld1
    Ld2
    Ld1sq
    Ld2sq
 end
%  
  methods
%    
    function obj = ACC_params()      
      obj.N=200; % number of fourier modes
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
      obj.Taabw = 0.1; % AABW export
      obj.rb=2e-4; % drag velocity
      obj.nu=1000; % Eddy viscosity
      obj.K=400; % Eddy diffusion
      obj.beta=1.5e-11; % coriolis parameter gradient
      obj.Ld1=sqrt(obj.rg*obj.H1)/abs(obj.f); % upper layer deformation radius
      obj.Ld2=sqrt(obj.rg*obj.H2)/abs(obj.f); % lower layer deformation radius
      obj.Ld1sq=1/obj.Ld1^2; % 1 over squared deformation radius for calculation
      obj.Ld2sq=1/obj.Ld2^2;
    end
%        
  end
end