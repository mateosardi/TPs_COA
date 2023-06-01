% Item 2 Observador sin integrador
% 
% Considerar que no puede medirse la corriente y sólo pueda medirse el ángulo, por lo
% que debe implementarse un observador. Obtener la simulación en las mismas condiciones que en
% el punto anterior, y superponer las gráficas para comparar.


clc;clear all;close all

% Definición de constantes
Laa=366e-6; J=5e-9; Ra=55.6; Bm=0; Ki=6.49e-3; Km=6.53e-3;
torque=1.15e-3;
t_etapa = 1e-5; tF = 0.6;
t=0:t_etapa:tF;
tiempo = round(tF/t_etapa);
ip(1)=0;

ent=(pi/2)*square(2*pi*t/0.2); %referencia
Tl=(torque/2)*square(2*pi*t/0.2)+torque/2; %torque aplicado periodicamente
ub=5*square(2*pi*t/0.6)+5; %accion de control de referencia

% Definición de matrices
% x1 = i, x2 = w, x3 = θ

A=[-Ra/Laa -Km/Laa 0; Ki/J -Bm/J 0; 0 1 0];
B=[1/Laa; 0; 0];
C=[0 0 1];
D=0;

% Condidiones iniciales
x0 = [0 0 0];
ia(1)=0;
w(1)=0;
theta(1)=0;
u(1) = 0;

% LQR
Q=diag([1 1/10000 1/40]);    R=0.1;

% Ganancia de realimentación
K = lqr(A, B, Q, R);

for i=1:1:tiempo
 x=[ia(i); w(i);theta(i)];
u(i)=-K*x; color='*-b'; 

ia_p=-Ra/Laa*ia(i)-Km/Laa*w(i)+1/Laa*u(i);
w_p=Ki/J*ia(i)-Bm/J*w(i)-Tl(i)/J;
theta_p=w(i);

ia(i+1)=ia(i)+t_etapa*ia_p;
w(i+1)=w(i)+t_etapa*w_p;
theta(i+1)=theta(i)+t_etapa*theta_p;

ip(i+1)=ia_p;

end


figure
plot(t,ia,'r');title('corriente');
xlabel('tiempo[s]');ylabel('angulo[rad]');legend('Ia')

figure
plot(ia,ip)
