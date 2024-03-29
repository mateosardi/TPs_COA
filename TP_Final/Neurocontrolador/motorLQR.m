clc;clear all;close all

% Definicion de constantes
LAA = 0.01; J = 0.00208;Ra = 1; B = 0.0011; LAF = 0.1238; RFF = 60; LFF = 60;
Vf = 220; % Tension de campo

torque=1.15e-3;
torque=0;
t_etapa = 1e-4; tF = 5;
t=0:t_etapa:tF;
tiempo = round(tF/t_etapa);
ip(1)=0;

% Definicion de matrices
A = [-Ra/LAA -LAF/LAA 0 0; -RFF/LFF 0 0 0; LAF/J 0 -B/J 0; 0 0 1 0];
B = [0; 1/LAA; 0; 0];
C = [0 0 0 1];
D = 0;

% Condiciones iniciales
% x = [ia(i); if(i); w(i); theta(i)]
x = [0; 0; 0; 10];
ia(1)=x(1);
i_f(1)=x(2);
w(1)=x(3);
theta(1)=x(4);
u(1) = 0;

% LQR
Q=diag([100 100 1/100000 1/100]);    R=10000;

% Ganancia de realimentacion
[K,Pr,e] = lqr(A, B, Q, R);

for i=1:1:tiempo
    u(i+1)=-K*x;
    x = mopdm2_motor(t_etapa,x,torque,u(i), Vf);
    ia(i+1) = x(1);
    i_f(i+1) = x(2);
    w(i+1) = x(3);
    theta(i+1) = x(4);
end

Pr
figure
subplot(2,2,1);
plot(t,u,'r');title('Accion de control');
xlabel('tiempo[s]');ylabel('Amplitud');
subplot(2,2,2);
plot(t,ia,'r');title('Corriente de armadura');
xlabel('tiempo[s]');ylabel('Corriente[A]');
subplot(2,2,3);
plot(t,theta,'r');title('Angulo \theta');
xlabel('tiempo[s]');ylabel('Angulo[rad]');
subplot(2,2,4);
plot(t,w,'r');title('Velocidad \omega');
xlabel('tiempo[s]');ylabel('Velocidad[rad/s]');


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
