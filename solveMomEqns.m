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
  Hb=params.Hb; % ridge height
  Wb=params.Wb; % ridge width
  Xb=params.Xb; % ridge longitude
  tau = params.tau;
  K = params.K;
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

  % generate Fourier coefficients for standing wave equations
  [c1,c2,c3,c4,c5,z1,z2] = calc_sw_coeffs(U,k,params);

  % standing wave solution
  psihat = zeros(2,N);
  psihat(2,:) = z2.*etabhat;
  psihat(1,:) = z1.*psihat(2,:);
  psi = N*real(ifft(psihat,[],2));

  %%% Check solution has converged
  tfs = -f*real(sum(1i*k.*psihat(2,:).*conj(etabhat),2));
  ifs = H2*Ld2sq*real(sum(1i*k.*psihat(2,:).*conj(psihat(1,:)),2));
  res=[tau/rho0 - rb*U(2) - tfs, ...
          -rb*U(2)+H2*K*Ld2sq*(U(1)-U(2))+ifs-tfs];
  if (sum(abs(res))>1e-10)
    error('Solution not converged');
  end
             
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
