%%%% Calculate Streamfunction and save to a file %%%%


%%% Initialize Arrays %%%%
y0 = yy_h;
Nx = 200;
Ny = length(y0);
U_all = zeros(2,Ny);
psi_all = zeros(2,Ny,Nx);
Psi = zeros(2,Ny,Nx);

%%% Loop through all lat's %%%
for j = 1:Ny
    params = ACC_params();
    H1 = params.H1;
    H2 = params.H2;
    %params.tau = .1.*sin(pi.*dy.*j./Ly).^2;
    %tau = .1.*sin(pi.*dy.*j./Ly).^2;
    [U, psi] = solveMomEqns (j);
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
 

%%% Mean the two layers 
Psi = squeeze((Psi(1,:,:)+Psi(2,:,:))/2);
%%% Save Streamfunction to a file %%%
save("Streamfunction.mat",'Psi')





