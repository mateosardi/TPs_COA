% clc; clear all; close all;

CI=0;
longx = 4;
dd = 1;
ho = 80;
% dW1 = zeros(ho,longx+1); dW2 = zeros(dd,ho+1);
% W1 = rand(ho,longx+1); W2 = rand(dd,ho+1);
% W1a=W1-0.25; W2a=W2-0.25;


tiempo_etapa=.1;
etapas=100;
xx1=0:etapas-1;%Angulo crítico 6.1*.2 rad aprox.
sx=[0;0;CI;0];
t=xx1;

% ref = pi/2*square(2*pi/50*t);
% plot(ref);

sal=zeros(4,1);   y1=zeros(1,etapas);y2=zeros(1,etapas);
vmax=0.2;u_a=0; Mw=zeros(4,etapas-1); costo=[];costo(1)=0;consigna=[];consigna(1)=0;
color='.-k';sigma_p=.01;sigma_tita=.1;
% color='.-r';sigma_p=.01;sigma_tita=.05;
% color='.-b';sigma_p=.01;sigma_tita=.02;
for k=1:etapas-1
    
    
   VX(:,k)=sx;
   entrada = [sx(1); sx(2);sx(3);sx(4)]+0*diag([sigma_p 0 sigma_tita 0])*randn(4,1);
   %Medición ruidosa
   y1(k)=entrada(3); %angulo del brazo
   y2(k)=entrada(1); %posición del carro   
   X = [ entrada; 1];     % bias
   s1 = W1a * X;
   y11=[pmntanh(s1); 1];
   s2 = W2a * y11;
   consigna(k)=s2;
   sal=mopdm2(tiempo_etapa,sx,consigna(k));
   sx=sal(:,1);
   costo_p(k)=indice_g(sal,consigna(k));
   costo(k+1)=costo(k)+indice_g(sx,consigna(k));
   u_a=consigna(k);
end

NPD=costo(etapas);
y1(k+1)=sal(3,1);
y2(k+1)=sx(1,1); %angulo del brazo
figure(10);
uo=consigna; To=5*6/109;
subplot(4,1,1),plot(To*xx1,y1,color),title('Evolución con neurocontrolador'),ylabel('Ángulo'),grid on,hold on;
subplot(4,1,2),plot(To*xx1,y2,color),ylabel('Posición'),grid on,hold on;
subplot(4,1,3),plot(To*xx1,costo,color),ylabel('Costo'),grid on,hold on; %,axis([0 etapas 0 6.1]);
subplot(4,1,4),plot(To*xx1(1:length(uo)),uo,color),ylabel('Acción de control'),xlabel('Tiempo [seg.]'),grid on,hold on;
