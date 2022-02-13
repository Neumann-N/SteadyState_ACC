% plot of schematic set up

% Load in data
params = ACC_params();
Lx=params.Lx;
Ly=params.Ly;
N=params.N;
f=params.f;
rb=params.rb;
rg=params.rg;
load("schematic.mat",'x','etab','psi','etab','H1','H2','k','K','U','nu','psihat')

%create plot
figure(4); hold on;
plot(x,etab,'b')
plot(x,(psi(1,:)/H1)+H1,'red')
plot(x,((psi(2,:))/H2)+H2)

%labels
annotation("textarrow",[.3 .4],[.3 .3],'String','U1 ')
annotation("textarrow",[.3 .4],[.7 .7],'String','U2 ')
annotation("textarrow",[.52 .52],[.14 .38],'String','Hb ')
annotation("textarrow",[.7 .7],[.14 .52],'String','H1 ')
annotation("textarrow",[.7 .7],[.55 .79],'String','H2 ')
%annotation("textarrow",[.4 .5],[.1 .1],'String','Wb ')
xlabel('x (m)')

legend('etab','psi1','psi2')