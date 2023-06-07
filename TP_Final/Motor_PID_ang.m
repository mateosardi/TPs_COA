% Implementación de un controlador PID en tiempo discreto para que el ángulo del motor permanezca en
% una referencia deseada.

clear;close all; clc;

% Variables de estado
% X = [ia, i_f, omega, theta]
X=[0; 0; 0; 0]; %valores iniciales

t_etapa=1e-3;
Ts=t_etapa;

titaRef=1; % Referencia deseada
tF=2; % Tiempo de simulación
Tl=0; % Torque aplicado
Vf = 12; % Tensión de campo

% Constantes del PID
Kp=80; Ki=300; Kd=10;

% Variables del PID
A1=((2*Kp*Ts)+(Ki*(Ts^2))+(2*Kd))/(2*Ts);
B1=(-2*Kp*Ts+Ki*(Ts^2)-4*Kd)/(2*Ts);
C1=Kd/Ts;

e=zeros(tF/t_etapa,1); % Vector de error
u=0;

ii=0;
for t=0:t_etapa:tF
    if t==0.5
        Tl=0
    end
    ii=ii+1;k=ii+2;
    X = modmotor_PID(t_etapa, X, u, Tl); % Función del modelo del motor
    e(k)=titaRef-X(4); %ERROR
    u = u + A1*e(k) + B1*e(k-1) + C1*e(k-2); % Acción de control utilizando PID
    
    % Saturación a 12V
    if u>12
        u = 12;
    end
    if u<-12
        u = -12;
    end
    
    i_a(ii)=X(1);% ia
    i_f(ii)=X(2);% if
    omega(ii)=X(3);% omega
    theta(k)=X(4);% theta
    acc(k)=u;
end

color_='r';
t=0:t_etapa:tF+2*t_etapa;
subplot(3,1,1);hold on;
plot(t,e,color_);title('Error');
subplot(3,1,2);hold on;
plot(t,acc,color_);title('Acción de control');
subplot(3,1,3);hold on;
plot(t,theta,color_);title('Ángulo \theta');
xlabel('Tiempo [Seg.]');
