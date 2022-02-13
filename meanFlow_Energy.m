function [Baroclinic_Conversion, Topo_Conversion, Creation_by_Wind, Destruction_by_Drag, Production_by_Buoyancy, Transient_eddy] = meanFlow_Energy(psi,U,tau,tfs,f,rg,rb,k,K,Taabw,rho0,psihat,N,Lx)    
    dpsihat(:,:)=1j*k.*psihat(:,:);
    Baroclinic_Conversion=-(U(1)-U(2)).*(f^2/rg)*real(sum(dpsihat(2,:).*conj(psihat(1,:)),2));
    Topo_Conversion=-U(2).*tfs;
    Creation_by_Wind=tau.*U(1)./rho0;
    Destruction_by_Drag=-rb.*U(2).^2;
    Production_by_Buoyancy=f.*Taabw.*(U(2)-U(1))./Lx;
    Transient_eddy=-K.*f.^2.*((U(1)-U(2)).^2)./rg;
end