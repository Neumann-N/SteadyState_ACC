%%%% Calculate Streamfunction and save to a file %%%%

clear all;

%%% Max wind stress
tau_max = 0.1;
K_max = 1000;
nu_max = 2000;

%%% Grid size
Nx = 256;
Ny = 128;

%%% Initialize model parameters
params = ACC_params();
H1 = params.H1;
H2 = params.H2;
Ly=params.Ly; % meridional domain size
Lx=params.Lx; % meridional domain size

%%% Modify parameters as needed
params.N = Nx;
params.Wb = 300000;
params.rb = 4e-4;
params.Hb = 1000;
params.rg = 1e-2;

%%% Initialize Arrays %%%%
dy = Ly/Ny;
dx = Lx/Nx;
yy = 0.5*dy:dy:Ly-0.5*dy;
xx = 0:dx:Lx-dx;
[XX,YY]=  meshgrid(xx,yy);

%%% Storage
U_all = zeros(2,Ny);
psi_all = zeros(2,Ny,Nx);
Psi = zeros(2,Ny,Nx);

%%% Loop through all lat's %%%
for j = 1:Ny

  params.tau = tau_max.*sin(pi.*yy(j)./Ly).^2;
  params.K = K_max*(params.tau/tau_max).^(1);
  params.nu = nu_max; %*params.tau/tau_max;
  [U, psi] = solveMomEqns (params);
  U_all(:,j) = U(:);
  psi_all(:,j,:) = psi(:,:);
  
end

%%% Calculate Mean Streamfunction %%%
for z = 1:Nx
  Psi(:,1,z) = -0.5*dy.*U_all(:,1);
end

for j = 2:Ny
  Psi(:,j,:) = Psi(:,j-1,:)-(0.5.*dy.*U_all(:,j-1))-(0.5.*dy.*U_all(:,j));  
end

%%% Total Streamfunction 
Psi(1,:,:) = (Psi(1,:,:) + psi_all(1,:,:)).*(H1);
Psi(2,:,:) = (Psi(2,:,:) + psi_all(2,:,:)).*(H2);
 


%%% Save Streamfunction to a file %%%
save("Streamfunction.mat",'Psi')





