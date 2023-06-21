clear; close all; clc;

% Variables de estado
% X = [ia, i_f, omega, theta]
X = [0; 0; 0; 0]; % Valores iniciales

t_etapa = 1e-3;
Ts = t_etapa;

titaRef = 1; % Referencia deseada
tF = 10; % Tiempo de simulación
Vf = 220; % Tensión de campo
Tl = 0.01; % Torque de carga

% Constantes del PID
Kp = 0.01; Ki = 0; Kd = 0;

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

% Dibujar el rotor y el estator utilizando patch
radio_rotor = 0.4; % Radio del rotor
radio_estator = 0.5; % Radio del estator

% Coordenadas del rotor y el estator (ejemplo)
theta_rotor = linspace(0, 2*pi, 100);
x_rotor = radio_rotor * cos(theta_rotor);
y_rotor = radio_rotor * sin(theta_rotor);
vertices_rotor = [x_rotor' y_rotor'];

theta_estator = linspace(0, 2*pi, 100);
x_estator = radio_estator * cos(theta_estator);
y_estator = radio_estator * sin(theta_estator);
vertices_estator = [x_estator' y_estator'];

% Definir las caras del rotor y el estator (ejemplo)
caras_rotor = 1:length(x_rotor);
caras_estator = 1:length(x_estator);

% Dibujar el rotor y el estator en la posición inicial
patch('Vertices', vertices_rotor, 'Faces', caras_rotor, 'FaceColor', 'r');
patch('Vertices', vertices_estator, 'Faces', caras_estator, 'FaceColor', 'b');

% Actualizar la posición del eje del motor en cada iteración
for i = 1:length(theta)
    posicion_eje = theta(i);
    
    cla; % Borra el gráfico anterior
    % Dibujar el rotor y el estator en la nueva posición del eje
    patch('Vertices', vertices_rotor + [posicion_eje, 0], 'Faces', caras_rotor, 'FaceColor', 'r');
    patch('Vertices', vertices_estator + [posicion_eje, 0], 'Faces', caras_estator, 'FaceColor', 'b');
    
    pause(0.1); % Pausa de 0.1 segundos para visualizar la posición actualizada
end
