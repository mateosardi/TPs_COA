clc;clear all;close all

% Definición de constantes
LAA = 0.01; J = 0.00208;Ra = 1; B = 0.0011; LAF = 0.1238; RFF = 60; LFF = 60;

torque=1.15e-3;
torque=0;
t_etapa = 1e-4; tF = 0.6;
t=0:t_etapa:tF;
tiempo = round(tF/t_etapa);
ip(1)=0;

% Definición de matrices
A=[-Ra/Laa -Km/Laa 0; Ki/J -Bm/J 0; 0 1 0];
B=[1/Laa; 0; 0];
C=[0 0 1];
D=0;

A =  [-Ra/LAA -LAF/LAA 0 0
      -RFF/LFF 0 0 0 
      LAF/J 0 -B/J 0 
      0 0 1 0 ];

B =  [0 0 
      1/LFF 0 
      0 0 
      0 0 ];



% Condidiones iniciales
% x = [ia(i); w(i); theta(i)]
x = [0; 0; 2];
ia(1)=x(1);
w(1)=x(2);
theta(1)=x(3);
u(1) = 0;

% LQR
Q=diag([100 1/900 1/4]);    R=1;

% Ganancia de realimentación
[K,Pr,e] = lqr(A, B, Q, R);

for i=1:1:tiempo
u(i+1)=-K*x;
x = mopdm2_motor(t_etapa,x,u(i));
ia(i+1) = x(1);
w(i+1) = x(2);
theta(i+1) = x(3);
end

Pr
figure
plot(t,u,'r');title('Accion de control');
xlabel('tiempo[s]');ylabel('Amplitud');legend('Ia')
figure
plot(t,ia,'r');title('corriente');
xlabel('tiempo[s]');ylabel('angulo[rad]');legend('Ia')
figure
plot(t,theta,'r');title('Ángulo \theta');
xlabel('tiempo[s]');ylabel('angulo[rad]');legend('Ia')
figure
plot(t,w,'r');title('Velocidad \omega');
xlabel('tiempo[s]');ylabel('angulo[rad]');legend('Ia')


% % Algoritmo aprendizaje Q
% Max_it=6; % Máximas iteraciones
% TM=100; % Cantidad de estados
% Mmax=3; % Valores que tiene uf
% t_etapa = 1e-4;
% color='.-k';
% tic; du=Mmax;
% etapas=50;umin=-10; umax=10;
% %%%Carga de datos
% rand('state',0);
% randn('state',0);
