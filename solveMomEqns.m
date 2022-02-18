%%%
%%% solveMomEqns.m
%%%
%%% Script to test ACC_optimizer: solves the mean flow/standing wave
%%% equations and evaluates the solution.
%%%

%%% Create parameter object with default parameters
function [U, psi] = solveMomEqns (tau)
params = ACC_params();

%%% Modify parameters as needed
%tau
%for q= [0.005 0.010 0.017 0.03 0.05 0.10 0.170 0.30]
%[0.001 0.0017 0.003 0.005 0.01 0.017 0.03 0.05 0.1 0.17 0.3 0.5]
%K
%for j= [25 50 100 146 200 250 300 310 320 395 405 430 445 495 505 535 540]
%[20 100 400 440 480 520 560 600 640 740]
%nu
%.*sin(pi.*y0./Ly).^2
%for y0 = yy_h
%for wsize= [40000.0 50000.0 60000.0 70000.0 80000.0 100000.0]
%for hsize= [400 450 500 600 650 750 1000 2000]
%for drr = [1e-4 1.5e-4 2e-4 2.5e-4 3e-4]
%Ly=params.Ly;
%params.tau = 0.1;

params.K = 600;
params.nu = 2000;
params.Wb=500000;
%params.Taabw = 0;

%%% Solve for mean flow
U = ACC_optimizer (params);

% extract model parameters
N=params.N; % number of fourier modes
Lx=params.Lx; % zonal domain size
Ly=params.Ly; % meridional domain size
H1=params.H1; % mean upper layer thickness
H2=params.H2; % mean bottom layer thickness
Hb=params.Hb; % ridge height
Wb=params.Wb; % ridge width
Xb=params.Xb; % ridge longitude
f=params.f; % coriolis parameter
%tau=params.tau; % wind stress
rho0=params.rho0; % density
rb=params.rb; % drag velocity    
g=params.g; % gravitational acceleration
rg=params.rg; % reduced gravity
nu=params.nu; % Eddy viscosity
K=params.K; % Eddy diffusion
beta=params.beta; % coriolis parameter gradient   
Ld1sq=params.Ld1sq; % 1 over squared deformation radius for calculation
Ld2sq=params.Ld2sq;

% grids in real and spectral space    
[x,k,etab,etabhat] = gen_grids (N,Lx,Hb,Xb,Wb);

% generate Fourier coefficients for standing wave equations
[c1,c2,c3,c4,c5,z1,z2] = calc_sw_coeffs(U,k,params);

% standing wave solution
psihat = zeros(2,N);
psihat(2,:) = z2.*etabhat;
psihat(1,:) = z1.*psihat(2,:);
psi = N*real(ifft(psihat,[],2));

%%% Compute IFS and TFS directly from psi1, psi2 and etab
%tfs = -f*real(sum(1i*k.*psihat(2,:).*conj(etabhat),2));
%ifs = H2*Ld2sq*real(sum(1i*k.*psihat(2,:).*conj(psihat(1,:)),2));

%%% Print a few things to check the solution
%tau/rho0;
%ifs;
%tfs;
%efs = H2*K*Ld2sq*(U(1)-U(2));
%drag = rb*U(2);
%tau/rho0-rb*U(2)-tfs;
%res=[tau/rho0 - rb*U(2) - tfs, ...
%         -rb*U(2)+H2*K*Ld2sq*(U(1)-U(2))+ifs-tfs];
             
%%% Plot the standng wave solution
%Tbc=H1.*(U(1)-U(2)).*Ly;
%Tbt=(H1+H2).*U(2).*Ly;
%Total=Tbc+Tbt;
%Tbt2=H2.*U(2).*Ly;
%figure(1)
%plot(x/1000,(f/g)*psi(1,:));
%ylabel('alpha');
%xlabel('wind stress');
%figure(i)
%plot(etab,x);
%ylabel('Isopycnal elevation (m)');
%xlabel('x (km)');

%save("Y_"+num2str(y0)+".mat",'Hb','Wb','tau','psihat','psi','etabhat','etab','H1','H2','k','K','ifs','tfs','nu','Tbc','Tbt','rb','res','U','Total','y0','x')
    %,"tfs","ifs","tau","K","nu","k","etabhat")
%end
%save("refernce_nu_"+num2str(tau)+"0.mat",'Hb','Wb','tau','psihat','psi','etabhat','etab','H1','H2','k','K','x','ifs','tfs','nu','Tbc','Tbt','rb','res','U')
%end
%end
%end
%end
%end