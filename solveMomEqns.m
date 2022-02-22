%%%
%%% solveMomEqns.m
%%%
%%% Script to test ACC_optimizer: solves the mean flow/standing wave
%%% equations and evaluates the solution.
%%%

%%% Create parameter object with default parameters
function [U, psi] = solveMomEqns (params)

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
