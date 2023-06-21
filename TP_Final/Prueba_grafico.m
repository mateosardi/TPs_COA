clear; close all; clc;

% Variables de estado
% X = [ia, i_f, omega, theta]
X = [0; 0; 0; 0]; % Valores iniciales

t_etapa = 1e-3;
Ts = t_etapa;

titaRef = 1; % Referencia deseada
tF = 4; % Tiempo de simulación
Vf = 220; % Tensión de campo
Tl = 0.01; % Torque de carga

% Constantes del PID
Kp = 1; Ki = 2; Kd = 0;

% Variables del PID
A1 = ((2 * Kp * Ts) + (Ki * (Ts^2)) + (2 * Kd)) / (2 * Ts);
B1 = (-2 * Kp * Ts + Ki * (Ts^2) - 4 * Kd) / (2 * Ts);
C1 = Kd / Ts;

e = zeros(tF / t_etapa, 1); % Vector de error
u = 0;

ii = 0;
for t = 0:t_etapa:tF
    ii = ii + 1; k = ii + 2;
    X = modmotor_PID(t_etapa, X, u, Vf, Tl); % Función del modelo del motor
    e(k) = titaRef - X(4); % ERROR
    u = u + A1 * e(k) + B1 * e(k - 1) + C1 * e(k - 2); % Acción de control utilizando PID
    
    % Saturación a 220 V
    if u > 220
        u = 220;
    end
    if u < -220
        u = -220;
    end
    
    i_a(ii) = X(1); % ia
    i_f(ii) = X(2); % if
    omega(ii) = X(3); % omega
    theta(ii) = X(4); % theta
    acc(ii) = u;
end

% Ajuste
i_a(ii + 1) = i_a(ii);
i_a(ii + 2) = i_a(ii);
i_f(ii + 1) = i_f(ii);
i_f(ii + 2) = i_f(ii);
theta(ii + 1) = theta(ii);
theta(ii + 2) = theta(ii);
omega(ii + 1) = omega(ii);
omega(ii + 2) = omega(ii);
acc(ii + 1) = acc(ii);
acc(ii + 2) = acc(ii);

% Crear figura y ejes
figure;
axes('XLim', [-1 1], 'YLim', [-1 1]);

% Dibujar el eje del motor
longitud_eje = 0.5; % Longitud del eje del motor
radio_rotor = 0.2; % Radio del rotor del motor
line([-longitud_eje/2, longitud_eje/2], [0, 0], 'Color', 'k', 'LineWidth', 2);

% Actualizar la posición del eje del motor en cada iteración
for i = 1:length(theta)
    posicion_eje = theta(i);
    
    cla; % Borra el gráfico anterior
    line([posicion_eje - longitud_eje/2, posicion_eje + longitud_eje/2], [0, 0], 'Color', 'k', 'LineWidth', 2); % Dibuja el nuevo eje del motor
    
    pause(0.1); % Pausa de 0.1 segundos para visualizar la posición actualizada
end
