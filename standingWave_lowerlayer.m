function [term1,term2,term3,term4,sum]= standingWave_lowerlayer(Ld1sq,psi,beta,nu,U,K,k)
psihat=fft(psi);
dpsihat=1j*k.*psihat;
d2psihat=-1.*k.^2.*psihat;
d3psihat=1j*-1*k.^3.*psihat;
d4psihat=k.^4.*psihat;
term1=-U(1).*(mean(d3psihat(1,:),2)+Ld1sq.*(mean(dpsihat(2,:),2)-mean(dpsihat(1,:),2)));
term2=mean(-dpsihat(1,:),2).*(beta+Ld1sq.*(U(1)-U(2)));
term3=(nu.*mean(d4psihat(1,:),2));
term4=(K.*Ld1sq.*(mean(d2psihat(2,:),2)-mean(d2psihat(1,:),2)));
sum=term1+term2+term3+term4;
end