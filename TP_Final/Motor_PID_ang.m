 % Implementación de un controlador PID en tiempo discreto para que el ángulo del motor permanezca en
% una referencia deseada.

% clear;close all; clc;

% Variables de estado
% X = [ia, i_f, omega, theta]
X=[0; 0; 0; CI]; %valores iniciales

t_etapa=1e-2;
Ts=t_etapa;

tF = 2; % Tiempo de simulación
Vf = 220; % Tensión de campo
% Tl = 0.01; % Torque de carga


% Variables del PID
A1=((2*Kp*Ts)+(Ki*(Ts^2))+(2*Kd))/(2*Ts);
B1=(-2*Kp*Ts+Ki*(Ts^2)-4*Kd)/(2*Ts);
C1=Kd/Ts;

e=zeros(tF/t_etapa,1); % Vector de error
u=0;

sigma_ia=.01;sigma_theta=.1;
ii=0;
for t=0:t_etapa:tF
    ii=ii+1;k=ii+2;
    X = modmotor_PID(t_etapa, X, u,Vf, Tl)+(0.5*diag([sigma_ia 0 0 sigma_theta])*randn(4,1))'; % Función del modelo del motor con ruido
    e(k)=titaRef-X(4); %ERROR
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
    omega(ii)=X(3);% omega
    theta(ii)=X(4);% theta
    acc(ii)=u;
end

% Ajuste
i_a(ii+1)=i_a(ii);
i_a(ii+2)=i_a(ii);
i_f(ii+1)=i_f(ii);
i_f(ii+2)=i_f(ii);
theta(ii+1) = theta(ii);
theta(ii+2) = theta(ii);
omega(ii+1)=omega(ii);
omega(ii+2)=omega(ii);
acc(ii+1)=acc(ii);
acc(ii+2)=acc(ii);

leyenda = sprintf('Torque (Tl = %.2f)', Tl);
% color_='r';
t=0:t_etapa:tF+2*t_etapa;
% subplot(3,1,1);hold on;
% plot(t,e);title('Error');
% % legend(leyenda);

% figure
subplot(2,2,1);
plot(t,acc);title('Accion de control'); hold on;
xlabel('tiempo[s]');ylabel('Amplitud [V]');
subplot(2,2,2);
plot(t,i_a);title('Corriente de armadura i_a'); hold on;
xlabel('tiempo[s]');ylabel('Corriente[A]');
subplot(2,2,3);
plot(t,omega);title('Velocidad \omega'); hold on;
xlabel('tiempo[s]');ylabel('Velocidad[rad/s]');
subplot(2,2,4);

plot(t,theta);title('Angulo \theta'); hold on;
% plot(t,mean(theta)*ones(size(t)),'-r'); hold on;
% plot(t,(mean(theta)+var(theta))*ones(size(t)),'r'); hold on;
% plot(t,(mean(theta)-var(theta))*ones(size(t)),'r'); hold on;
xlabel('tiempo[s]');ylabel('Angulo[rad]'); 


fprintf('Media: %.4f\n',mean(theta))
fprintf('Varianza: %.4f\n',var(theta))

