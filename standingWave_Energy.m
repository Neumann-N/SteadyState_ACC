function [Baroclinic_Conversion,Topo_Conversion,Drag_Dissipation,Viscos_Dissipation,Conversion_to_Transient_Eddy,totals] = standingWave_Energy(psi,U,tfs,f,rg,rb,nu,H1,H2,k,K,psihat,N)
    dpsihat(:,:)=1j.*k.*psihat(:,:);
%     d2psihat(:,:)=(1j.*k).^2.*psihat(:,:);
    Baroclinic_Conversion=(U(1)-U(2)).*(f^2/rg).*sum(real((dpsihat(2,:).*conj(psihat(1,:)))),2);
    Topo_Conversion=U(2).*tfs;
    Drag_Dissipation(:)=-rb.*sum(real(k.^2.*abs(psihat(2,:)).^2),2);
%     Drag_Dissipation(:)=-rb.*sum(real(dpsihat(2,:).*conj(dpsihat(2,:))),2);
    Viscos_Dissipation(:)=-nu.*(H1.*((sum(real(k.^4.*abs(psihat(1,:)).^2),2)))+(H2*sum(real(k.^4.*abs(psihat(2,:)).^2),2)));
%     Viscos_Dissipation(:)=-nu.*sum(real( H1*d2psihat(1,:).*conj(d2psihat(1,:)) + H2*d2psihat(2,:).*conj(d2psihat(2,:)) ),2);
    Conversion_to_Transient_Eddy(:)=-K.*f^2.*(((sum(real(k.^2.*abs((psihat(1,:)-psihat(2,:))).^2),2))))/rg;
%     Conversion_to_Transient_Eddy(:)=-K.*f^2.*(((sum(real((dpsihat(1,:)-dpsihat(2,:)).*conj(dpsihat(1,:)-dpsihat(2,:))),2))))/rg;
    totals=Baroclinic_Conversion+Topo_Conversion+Drag_Dissipation+Viscos_Dissipation+Conversion_to_Transient_Eddy;
end