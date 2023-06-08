% Implementación de un controlador PID en tiempo discreto para que la velocidad del motor permanezca en
% una referencia deseada.

clear;close all; clc;

% Variables de estado
% X = [ia, i_f, omega, theta]
X=[0; 0; 0; 0]; %valores iniciales

t_etapa=1e-3;
Ts=t_etapa;

omegaRef=1000*2*pi/60; % Referencia deseada en [rad/s]
tF=2; % Tiempo de simulación
Tl=0.001; % Torque aplicado
Vf = 220; % Tensión de campo

% Constantes del PID
Kp=15; Ki=50; Kd=0;

% Variables del PID
A1=((2*Kp*Ts)+(Ki*(Ts^2))+(2*Kd))/(2*Ts);
B1=(-2*Kp*Ts+Ki*(Ts^2)-4*Kd)/(2*Ts);
C1=Kd/Ts;

e=zeros(tF/t_etapa,1); % Vector de error
u=0;

ii=0;
for t=0:t_etapa:tF
    ii=ii+1;k=ii+2;
    X = modmotor_PID(t_etapa, X, u, Vf, Tl); % Función del modelo del motor
    e(k)=omegaRef-X(3); %ERROR
    u = u + A1*e(k) + B1*e(k-1) + C1*e(k-2); % Acción de control utilizando PID
    
    % Saturación a 220 V
    if u>220
        u = 220;
    end
    if u<-220
        u = -220;
    end
    
    i_a(ii)=X(1);% ia
    i_f(ii)=X(2);% if
    omega(k)=X(3);% omega
    theta(ii)=X(4);% theta
    acc(k)=u;
end

color_='r';
t=0:t_etapa:tF+2*t_etapa;
subplot(3,1,1);hold on;
plot(t,e,color_);title('Error');
subplot(3,1,2);hold on;
plot(t,acc,color_);title('Acción de control');
subplot(3,1,3);hold on;
plot(t,omega,color_);title('Velocidad \omega');
xlabel('Tiempo [Seg.]');
