clearvars
clc
params = ACC_params();
N=params.N; % number of fourier modes
H1=params.H1; % mean upper layer thickness
H2=params.H2; % mean bottom layer thickness
f=params.f; % coriolois parameter
rb=params.rb; % drag velocity    
rg=params.rg;
Lx=params.Lx;
beta=params.beta;
%Hb=params.Hb;
Xb=params.Xb;
%Wb=params.Wb;
Ld1sq=params.Ld1sq; % 1 over squared deformation radius for calculation
Ld2sq=params.Ld2sq;
savealpha=0;
%w = dir("nu1000_*_K400.mat");
w = dir("Wb40000_*.mat");
for j=1:21
filename2 = w(j).name;
load(filename2,'Hb','Wb','tau','psihat','psi','etabhat','etab','H1','H2','k','K','ifs','tfs','nu','Tbc','Tbt','rb','res')

%Heights=[500 1000 1500];
%width=[40000 80000 150000 300000 600000];

% Solving for Alpha 
EPE=.5.*(f.^2)./rg.*((sum(real(abs((psihat(1,:)-psihat(2,:))).^2),2)));
EKE=(.5.*H1.*sum(real(k.^2.*abs(psihat(1,:)).^2),2)) + (.5.*H2.*sum(real(k.^2.*abs(psihat(2,:)).^2),2));
%E2=EPE;
E2=EPE+EKE;
bouyancyf=(rg./(H1+H2)).^.5;
%alpha1=ifs./E1./abs(f).*bouyancyf.*(H1+H2);
alpha=ifs./E2./abs(f).*bouyancyf.*(H1+H2);
savealpha=[savealpha alpha];
end
%byhand=(2.*real(sum(1i*k.*psihat(2,:).*conj(psihat(1,:)),2))./((sum(real(abs((psihat(1,:)-psihat(2,:))).^2),2))))./abs(f);
tau= [0.01 0.017 0.03 0.05 0.1 0.17 0.3 0.5];
%savealpha(2:7:-1)
figure(6); hold on;
%P1=semilogx(tau,alpha1,"o","Color",'r');
P1=semilogx(tau,savealpha(1:8),'-o',"Color",'b');
P2=semilogx(tau,savealpha(8:15),'-o',"Color",'r');
P3=semilogx(tau,savealpha(15:22),'-o',"Color",'g');
set(gca,'xscale','log')


legend([P1 P2 P3],'Hb=500m',"Hb=1000m","Hb=1500m")
title('Wb = 600000m')
xlabel("tau (N/m^2)")
ylabel("alpha")


%term1=nu.*H2./U(2).*sum(real(k.^4.*abs(psihat(2,:)).^2),2);
%term2=K.*Ld1sq.*H1./U(1).*sum(real(k.^2.*abs(psihat(1,:)).^2),2);

%for p(1) equation 2
%term1=U(1).*Ld1sq.*real(sum(1i*k.*psihat(1,:).*conj(psihat(2,:)),2));
%term2=nu.*(sum(real(k.^4.*abs(psihat(1,:)).^2),2));
%term31=-K.*Ld1sq.*(sum(real(1i*k.*(psihat(1,:)).*conj(1i.*k.*psihat(2,:))),2))
%term32=-K.*Ld1sq.*(-(sum(real(k.^2.*abs(psihat(1,:)).^2),2)));

%for p(2) equation 3
%term1=f.*U(2).*real(sum(1i*k.*psihat(2,:).*conj(etabhat),2))./H2;
%term2=rb.*sum(real(k.^2.*abs(psihat(2,:)).^2),2)./H2;
%term32=nu.*(sum(real(k.^4.*abs(psihat(2,:)).^2),2));
%term31=nu.*(sum(real(k.^4.*abs(psihat(1,:)).^2),2).*U(2).*H1./(U(1).*H2));
%term4=-K.*Ld2sq.* ...
%    (sum(real( 1i*k.*(psihat(1,:).*conj(1i*k.*psihat(2,:))) ),2) ...
%    -sum(real( k.^2.*abs(psihat(2,:)).^2 ),2) ...
%    +U(2).*(sum(real( 1i*k.*(psihat(1,:).*conj(1i*k.*psihat(2,:))) ),2))./U(1) ...
%    -U(2).*sum(real( k.^2.*abs(psihat(1,:)).^2 ),2)./U(1));

%for (5)
%dx=Lx/N;
%x = -Lx/2:dx:Lx/2-1;
%etab = Hb * exp(-((x-Xb)/(.5*Wb)).^2);
%term1=sum(real(psi(2,:)),2);
%term2=(psi(2,:)),psi(1,:).*(U(2)-(beta.*Ld1sq))./U(1);
%term1=U(1)-U(2);

%figure(4); hold on;
%P1=semilogx(tau,alpha1,"o","Color",'r');
%P2=semilogx(tau,alpha2,"o","Color",'b');
%set(gca,'xscale','log')
%scatter(term1,term2,'b')
%scatter(psi(2,:),-psi(1,:).*(U(2)-(beta.*Ld1sq.^-1)./U(1)-(f.*etab(1,:).*Ld1sq.^-1)./H2./260000),'blue','filled')
%A=scatter(tau,psi(2,:),'green');
%B=scatter(tau,psi(1,:).*(U(2)-(rb.*Ld1sq.^(-1))./U(1)./H1./Lx),'magenta');
%scatter(psi(2,:),U(1)/U(2).*psi(1,:))
%P1=semilogx(tau,term1,"o","Color",'b');
%P2=semilogx(tau,term2,"o","Color",'r');
%P3=semilogx(tau,term31./term1,"o","Color",'y');
%P4=semilogx(tau,term32./term1,"o","Color",'g');
%P5=semilogx(tau,term4./term1,"o","Color",'k');
%[Baroclinic_Conversion,Topo_Conversion,Drag_Dissipation,Viscos_Dissipation,Conversion_to_Transient_Eddy] = standingWave_Energy(psi,U,tfs,f,rg,rb,nu,H1,H2,k,K,psihat,N)
%total=Baroclinic_Conversion+Topo_Conversion+Drag_Dissipation+Viscos_Dissipation+Conversion_to_Transient_Eddy
%end

%xlabel("v*H2*psi2xx^2/U2")
%ylabel("K*Ld1^-2*H1*psi1x^2/U1")
%legend([P2 P3 P4 P5],"rb*psi2x^2/H2/'U2*f/H2*etab*psi2x'",'nu*(psi2xx^2+U2*H1/(U1*H2)*psi1xx^2)/U2*f/H2*etab*psi2x','K*Ld2*(psi1x*psi2x-psi2x^2+U2/U1*psi1x*psi2x-U2/U1*psi1x^2)/U2*f/H2*etab*psi2x','total');
%legend([P1 P2],'potential energy only',"Both kenetic and potential")
%title('nu = ')

%total=term1+term2+term3+term4
%term1
%term2
%term31
%term4