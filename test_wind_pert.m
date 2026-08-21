%%%
%%% test_wind_pert.m
%%%
%%% Test whether the barotropic streamfunction from a wind perturbation
%%% can account for the change in IFS.
%%%

load('Streamfunction_tau0.1.mat');
Psi1 = psi_all;
PsiBT1 = squeeze((psi_all(1,:,:) + psi_all(2,:,:)) / 2);
PsiBC1 = squeeze((psi_all(1,:,:) - psi_all(2,:,:)) / 2);
load('Streamfunction_tau0.15.mat');
Psi2 = psi_all;
PsiBT2 = squeeze((psi_all(1,:,:) + psi_all(2,:,:)) / 2);
PsiBC2 = squeeze((psi_all(1,:,:) - psi_all(2,:,:)) / 2);

% load('Streamfunction.mat');

%%% Grid sizes
Ny = size(Psi,2);
Nx = size(Psi,3);

%%% Initialize Arrays %%%%
dy = Ly/Ny;
dx = Lx/Nx;
yy = 0.5*dy:dy:Ly-0.5*dy;
xx = 0:dx:Lx-dx;
[XX,YY]=  meshgrid(xx,yy);


% extract model parameters
N=params.N; % number of fourier modes
Lx=params.Lx; % zonal domain size
Hb=params.Hb; % ridge height
Wb=params.Wb; % ridge width
Xb=params.Xb; % ridge longitude
tau = params.tau;
Ky = params.Ky;
rho0 = params.rho0;
rb = params.rb;
H2 = params.H2;
H1 = params.H1;
rg = params.rg;
f = params.f;
Ld1=sqrt(rg*H1)/abs(f); %%% upper layer deformation radius
Ld2=sqrt(rg*H2)/abs(f); %%% lower layer deformation radius
Ld1sq=1/Ld1^2; 
Ld2sq=1/Ld2^2;

% grids in real and spectral space    
[x,k,etab,etabhat] = gen_grids (N,Lx,Hb,Xb,Wb);

tfs = zeros(Ny,1);
ifs = zeros(Ny,1);
for j=1:Ny
  
  % standing wave solution  
  psi1hat = fft(squeeze(Psi1(:,j,:)),[],2)/Nx;  
  psi2hat = fft(squeeze(Psi2(:,j,:)),[],2)/Nx;  

  PsiBT1hat = fft(PsiBT1(j,:),[],2)/Nx;  
  PsiBC1hat = fft(PsiBC1(j,:),[],2)/Nx;  
  PsiBT2hat = fft(PsiBT2(j,:),[],2)/Nx;  
  PsiBC2hat = fft(PsiBC2(j,:),[],2)/Nx;  

  %%% TFS and IFS calculation using psi1 and psi2
%   tfs(j) = -f*real(sum(1i*k.*psi1hat(2,:).*conj(etabhat),2));
%   ifs(j) = H2*Ld2sq*real(sum(1i*k.*psi1hat(2,:).*conj(psi1hat(1,:)),2));
%   tfs(j) = -f*real(sum(1i*k.*psi2hat(2,:).*conj(etabhat),2));
%   ifs(j) = H2*Ld2sq*real(sum(1i*k.*psi2hat(2,:).*conj(psi2hat(1,:)),2));
%   tfs(j) = -f*real(sum(1i*k.*psi2hat(2,:).*conj(etabhat),2));
%   ifs(j) = H2*Ld2sq*real(sum(1i*k.*psi2hat(2,:).*conj(psi1hat(1,:)),2));

  %%% TFS and IFS calculation using psiBT and psiBC
%   tfs(j) = -f*real(sum(1i*k.*(PsiBT1hat-PsiBC1hat).*conj(etabhat),2));
%   ifs(j) = 2*H2*Ld2sq*real(sum(1i*k.*PsiBT1hat.*conj(PsiBC1hat),2));
%   tfs(j) = -f*real(sum(1i*k.*(PsiBT2hat-PsiBC2hat).*conj(etabhat),2));
%   ifs(j) = 2*H2*Ld2sq*real(sum(1i*k.*PsiBT2hat.*conj(PsiBC2hat),2));
  tfs(j) = -f*real(sum(1i*k.*(PsiBT2hat-PsiBC1hat).*conj(etabhat),2));
  ifs(j) = 2*H2*Ld2sq*real(sum(1i*k.*PsiBT2hat.*conj(PsiBC1hat),2));

end


figure(1);
clf
plot(yy,ifs);
hold on;
plot(yy,tfs);
plot(yy,-rb*U_all(2,:));
hold off;