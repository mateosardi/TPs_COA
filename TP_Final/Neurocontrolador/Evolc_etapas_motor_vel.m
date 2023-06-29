% clc; clear all; 
close all;
TamanioFuente=14;
CI=0;
titaRef=0;
color='.r';colorc='r';
tiempo_etapa =t_etapa;
Ts = t_etapa;

torque = 0.1;
etapas=1500;
xx1=0:etapas-1;
sx=[0;0;0;CI];
t=xx1;

sal=zeros(4,1);   y4=zeros(1,etapas);y1=zeros(1,etapas);
u_a=0; Mw=zeros(4,etapas-1); costo=[];costo(1)=0;consigna=[];consigna(1)=0;
color='.-k';

sigma_ia=.01;sigma_theta=.1;sigma_omega=.1;

u=0;
for k=1:etapas-1
    VX(:,k)=sx;
    entrada = [sx(1); sx(2);sx(3);sx(4)]+15*diag([sigma_ia 0 sigma_omega sigma_theta])*randn(4,1);
    %Medicion ruidosa
    y1(k)=entrada(1); % corriente armadura
    y2(k)=entrada(2); % corriente de campo
    y3(k)=entrada(3); % velocidad angular
    y4(k)=entrada(4); % angulo del motor
    X = [entrada.*V_max; 1]; % bias
    s1 = W1a * X;
    y11=[pmntanh(s1); 1];
    s2 = W2a * y11;
    
    consigna(k)=s2;

    sal=mopdm2_motor(tiempo_etapa,sx,torque,consigna(k),Vf);
    wp(k) = sal(5);
    sal = sal(1:4);
    sx=sal(:,1);
    costo_p(k)=indice_g_motor(sal,consigna(k));
    costo(k+1)=costo(k)+indice_g_motor(sx,consigna(k));
    u_a=consigna(k);    
end


NPD=costo(etapas);

y1(k+1)=sal(1,1); % corriente armadura
y2(k+1)=sal(2,1); % corriente de campo
y3(k+1)=sal(3,1); % velocidad angular
y4(k+1)=sal(4,1); % angulo motor
wp(k+1) = wp(k);


uo=consigna; To=t_etapa;
figure;
plot(To*xx1(1:length(uo)),uo,color),ylabel('Accion de control'),xlabel('Tiempo [seg.]'),grid on,title('Accion de control');

figure;
subplot(2,2,1);
plot(To*xx1(1:length(uo)),uo);
ylabel('Accion de control'),xlabel('Tiempo [seg.]'),title('Accion de control');
subplot(2,2,2);
plot(To*xx1,y1);title('Corriente de armadura ia'); hold on;
xlabel('tiempo[s]');ylabel('Corriente[A]');
subplot(2,2,3);
plot(To*xx1,y3);title('Velocidad \omega'); hold on;
xlabel('tiempo[s]');ylabel('Velocidad[rad/s]'); 
subplot(2,2,4);
plot(To*xx1,y4);title('Angulo \theta'); hold on;
xlabel('tiempo[s]');ylabel('Angulo[rad]');

t_theta = [To*xx1;y4]';
t_omega = [To*xx1;y3]';
t_wp = [To*xx1;wp]';
