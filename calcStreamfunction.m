%%%% Calculate Streamfunction and save to a file %%%%


%%% Initialize Arrays %%%%
y0=yy_h;
Nx=200;
Ny = length(y0);
U_all = zeros(2,Ny);
psi_all = zeros(2,Ny,Nx);
Psi = zeros(2,Ny,Nx);
delY = y0(2)-y0(1);

%%% Loop through all lat's %%%
for j = 1:Ny
    params = ACC_params();
    H1=params.H1;
    params.tau = .1.*sin(pi.*delY.*Ny./Ly).^2;
    tau=params.tau;
    [U, psi] = solveMomEqns (tau);
    U_all(:,j) = U;
    psi_all(:,j,:) = psi;
end
U_all=squeeze(U_all(1,:)-U_all(2,:))*H1;
%%% Calculate Streamfunction %%%
Psi(:,1,:) = -0.5*delY*U_all(1).*ones(2,200);
for j=2:Ny
  
  %%% Mean streamfunction
  Psi(:,j,:) = Psi(:,j-1,:)-0.5*delY*U_all(j-1)-0.5*delY*U_all(j);
  %%% Standing wave stream function 
  Psi(:,j,:) = Psi(:,j,:) + psi_all(:,j,:);
end

%%% Mean the two layers 
Psi=squeeze((Psi(1,:,:)+Psi(2,:,:))/2);

%%% Save Streamfunction to a file %%%
save("Streamfunction.mat",'Psi')

