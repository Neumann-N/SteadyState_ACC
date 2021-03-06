function [term1,term2,term3,term4,term5,sum]= standingWave_upperlayer(Ld2sq,psi,beta,nu,U,f,H2,K,k,rb,etabhat)
psihat=fft(psi);
dpsihat=1j*k.*psihat;
d2psihat=-1.*k.^2.*psihat;
d3psihat=1j*-1*k.^3.*psihat;
d4psihat=k.^4.*psihat;
term1=-U(2).*(mean(d3psihat(2,:),2)+Ld2sq.*(mean(dpsihat(1,:),2)-mean(dpsihat(2,:),2))+(f.*mean(etabhat)./H2));
term2=mean(-dpsihat(2,:),2).*(beta+Ld2sq.*(U(2)-U(1)));
term3=(nu.*mean(d4psihat(2,:),2));
term4=(K.*Ld2sq.*(mean(d2psihat(1,:),2)-mean(d2psihat(2,:),2)));
term5=-rb*mean(d2psihat(2,:),2)/H2;
sum=term1+term2+term3+term4+term5;
end
