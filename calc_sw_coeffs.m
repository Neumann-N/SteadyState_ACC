%%%
%%% calc_sw_coeffs.m
%%%
%%% Compute complex coeffients from standing wave equations in spectral
%%% space.
%%%
function [c1,c2,c3,c4,c5,z1,z2] = calc_sw_coeffs(U,k,params)

  % extract model parameters    
  H2=params.H2; % mean bottom layer thickness  
  f=params.f; % coriolois parameter  
  rb=params.rb; % drag velocity    
  nu=params.nu; % Eddy viscosity
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
  
  % compute coefficients
  c1=-U(1)*k.^2+beta-Ld1sq*U(2)+1i*k.^3*nu+1i*k*Ld1sq*K;
  c2=U(1)*Ld1sq-1i*k*Ld1sq*K;
  c3=U(2)*Ld2sq-1i*k*Ld2sq*K;
  c4=-U(2)*k.^2+beta-Ld2sq*U(1)+(1i*k*rb)/H2+1i*k.^3*nu+1i*k*Ld2sq*K;
  c5=(U(2)*f)/H2;
  z1=-c2./c1;
  z2=-c5./(c3.*z1+c4);

end

