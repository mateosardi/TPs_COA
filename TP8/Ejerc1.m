% Programación dinámica.
% Apredizaje Q.
% Para el Ejemplo 7.3.1.
% x(k+1)=x(k)+(2-2*x(k)+5/4*x(k)^2-1/4*x(k)^3)*u(k);
% con el funcional de costo J=sum((2+u(k))*exp(-x(k)));
%Autor JAP
clear all;clc;close all;
TM=200; Mmax=16; color='.-k';

tic; du=Mmax; etapas=6; xmin=0; xmax=3;umin=-1; umax=1;
%%%Carga de datos
rand('state',0); equis1=3*(rand(TM,1)); tiempo=ceil((etapas-1)*rand(TM,1));
M_est = [tiempo,equis1]';
Au=(umax-umin)/(Mmax-1);
for i=1:Mmax
uf(i)=umin+Au*(i-1);
end
CI=2;vfinal=1; Q = zeros(TM,du);
J=zeros(1,TM); sal(1)=CI;costo(1)=0;
Ya=zeros(1,TM);
for k=1:etapas-1
entrada = [k; sal(k)];
consigna(k) = pol_tab_mu1(entrada,M_est,Ya);
sal(k+1)=mopdm(k,sal(k),consigna(k));
costo(k+1)=costo(k)+indice(k,sal(k),consigna(k));
end
costo(k+1)=costo(k+1)+abs(sal(k+1)-vfinal);
evoluc(1)=costo(etapas); costo(etapas)
m=zeros(TM,du);
for iterac=1:10
for iq=1:TM
k=M_est(1,iq);
x=M_est(2,iq);
%Recorre todas las acciones
for acc=1:du
xy=mopdm(k,x,uf(acc));
m(iq,acc)=m(iq,acc)+1;
gama=10/(10+m(iq,acc));
if k<etapas-1
[aux lugar]=min(abs(((k+1)-M_est(1,:))/max(M_est(1,:)))+abs((xy-M_est(2,:))/max(M_est(2,:))));
Q(iq,acc)=(1-gama)*Q(iq,acc)+gama*(indice(k,x,uf(acc))+J(lugar));
else
Q(iq,acc)=(1-gama)*Q(iq,acc)+gama*(indice(k,x,uf(acc))+abs(xy-vfinal));
end
end
end
for iq=1:TM
[val lugar]=min(Q(iq,:));
J(iq)=val;
Ya(iq)=uf(lugar);
end
sal(1)=CI;
costo=0;
for k=1:etapas-1
entrada = [k; sal(k)];
consigna(k) = pol_tab_mu1(entrada,M_est,Ya);
sal(k+1)=mopdm(k,sal(k),consigna(k));
costo(k+1)=costo(k)+indice(k,sal(k),consigna(k));
end
costo(k+1)=costo(k+1)+abs(sal(k+1)-vfinal);
evoluc(iterac)=costo(etapas); costo(etapas)
 plot(sal),pause
end
toc
