%%%% Calculate Streamfunction and save to a file %%%%

% clear all;

%%% Max wind stress
tau_max = 0.3;
K_max = 400;
nu_max = 2000;

%%% Grid size
Nx = 256;
Ny = 128;
m1km = 1000;

%%% Initialize model parameters
params = ACC_params();
H1 = params.H1;
H2 = params.H2;
Ly=params.Ly; % meridional domain size
Lx=params.Lx; % meridional domain size

%%% Modify parameters as needed
params.N = Nx;
params.Wb = 150000;
params.rb = 4e-4;
params.Hb = 1000;
params.Ky = 80;
params.Hrms = 0;
params.nu1 = 4000;
params.nu2 = 4000;
params.Xb = 1000000;
% params.Xb = Lx/4;
% params.Xb2 = 3*Lx/4;
params.double_ridge = false;

%%% Initialize Arrays %%%%
dy = Ly/Ny;
dx = Lx/Nx;
yy = 0.5*dy:dy:Ly-0.5*dy;
xx = 0:dx:Lx-dx;
[XX,YY]=  meshgrid(xx,yy);

%%% Generate bathymetry
rng(6); %%% Fix random seed
etab = params.Hb*exp(-((XX-params.Xb)/params.Wb).^2);
if (params.double_ridge)
    etab = etab+params.Hb*exp(-((XX-params.Xb2)/params.Wb).^2);
end
etab = etab + genBathy(400000,0,params.Hrms,-3,Nx,Ny,Lx,Ly)';

%%% Storage
U_all = zeros(2,Ny);
psi_all = zeros(2,Ny,Nx);
Psi = zeros(2,Ny,Nx);

%%% Loop through all lat's %%%
for j = 1:Ny

  params.tau = tau_max.*sin(pi.*yy(j)./Ly).^2;
  params.K = K_max;%*(params.tau/tau_max).^(1);
  params.nu1 = nu_max; %*params.tau/tau_max;
  params.nu2 = nu_max; %*params.tau/tau_max;  
  params.etab = etab(j,:);
  
  [U, psi] = solveMomEqns (params);
  U_all(:,j) = U(:);
  psi_all(:,j,:) = psi(:,:);

  % [x_tmp,k_tmp,etab_tmp,etabhat_tmp] = gen_grids (Nx,Lx,etab);
  % etab(:,j) = etab_tmp;

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

Ttot = sum(H1*U_all(1,:)+H2*U_all(2,:))*dy
Tbt = sum(U_all(2,:)*(H1+H2))*dy
Tbc = Ttot - Tbt


figure(2);
clf;
axes;
set(gca,'Position',[0.07 0.1 0.91 0.83])
pcolor(XX/m1km,YY/m1km,-params.H1 + squeeze(psi_all(2,:,:)-psi_all(1,:,:))*params.f/params.rg);
shading interp;
hold on;
[C,h]=contour(XX/m1km,YY/m1km,params.H1+params.H2-etab,[500 1500 2500 3250 3500 3900],'EdgeColor',[.5 .5 .5],'LineWidth',1);
hold off;
xlabel('Longitude $x$ (km)','interpreter','latex');
ylabel('Latitude $y$ (km)','interpreter','latex');
colorbar;
box on;
caxis([-2500 -1500]);
colormap(gca,haxby(22));
set(gca,'FontSize',14);
title('Mean sea surface height (m)');


 
figure(1);
clf;
set(gcf,'Position',[400 400 800 800]);

subplot('Position',[0.09 0.55 0.86 0.42])
pcolor(XX/m1km,YY/m1km,squeeze(Psi(1,:,:))/H1*params.f/params.g);
shading interp;
hold on;
[C,h]=contour(XX/m1km,YY/m1km,params.H1+params.H2-etab,[500 1500 2500 3250 3500 3900],'EdgeColor',[.5 .5 .5],'LineWidth',1);
hold off;
% xlabel('Longitude $x$ (km)','interpreter','latex');
ylabel('Latitude $y$ (km)','interpreter','latex');
handle = colorbar;
set(handle,'Position',[0.95 0.55 0.015 0.42]);
set(gca,'Position',[0.09 0.55 0.84 0.42])
box on;
caxis([0 2]);
colormap(gca,haxby(22));
set(gca,'FontSize',14);
title('Mean sea surface height (m)');
annotation('textbox',[0 0.5 0.05 0.03],'String','(a)','interpreter','latex','FontSize',16,'LineStyle','None');


load('ModelRefDiags.mat','ssh','eta','XX_h');
xx_h = XX_h(:,1);
ssh1D = ssh(:,end/2);
eta1D = eta(:,end/2);
% ssh1D = mean(ssh,2);
% eta1D = mean(eta,2);
Psi1D = squeeze(psi_all(:,Ny/2,:));
% Psi1D = squeeze(mean(psi_all,2));
subplot('Position',[0.09 0.06 0.84 0.42])
% plot(xx/m1km,squeeze(Psi(:,Ny/2,:)-repmat(mean(Psi(:,Ny/2,:),3),[1 1 Nx])));
plot(xx/m1km,1000*Psi1D(1,:)*params.f/params.g,'Color',[48 128 238]/256,'LineWidth',2);
hold on;
plot(xx_h/m1km,1000*ssh1D,'Color',[48 128 238]/256,'LineWidth',2,'LineStyle','--');
plot(xx/m1km,-params.H1 + (Psi1D(2,:)-Psi1D(1,:))*params.f/params.rg,'Color',[24 60 139]/256,'LineWidth',2);
plot(xx_h/m1km,eta1D,'Color',[24 60 139]/256,'LineWidth',2,'LineStyle','--');
plot(xx/m1km,-H1-H2+etab(Ny/2,:),'Color',[139 69 19]/256,'LineWidth',2);
hold off
set(gca,'XLim',[0 Lx/m1km]);
xlabel('Longitude $x$ (km)','interpreter','latex');
ylabel('Elevation $z$ (m)','interpreter','latex');
box on;
set(gca,'FontSize',14);
legend('Sea surface height \times 1000 (theory)','Sea surface height \times 1000 (model)','Isopycnal elevation (theory)','Isopycnal elevation (model)','Bathymetry','Location','SouthEast');
annotation('textbox',[0 0 0.05 0.03],'String','(b)','interpreter','latex','FontSize',16,'LineStyle','None');

%%% Save Streamfunction to a file %%%
save("Streamfunction.mat",'Psi','params','U_all','psi_all')





