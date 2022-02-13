%%%
%%% gen_grids.m
%%%
%%% Generates coordinates and bathymetry in real and spectral space.
%%%
function [x,k,etab,etabhat] = gen_grids (N,Lx,Hb,Xb,Wb)
%​
  % grids in real and spectral space    
  dx=Lx/N;
  x = -Lx/2:dx:Lx/2-1;
  nn = zeros(1,N);
  nn(1:N/2) = 0:1:N/2-1;
  nn(N/2+1:N) = -N/2:1:-1;
  k = 2*pi*nn/Lx; 
  %n=5;
%​
  % bathymetry 
  etab=Hb.*exp(-(x./Wb).^2); % calculates ridge normal
  %etab = Hb * cos(2.*pi.*n.*x./Lx); %for ridge effects data
  etabhat=(fft(etab))/N; % put etab in spectral space
  
end