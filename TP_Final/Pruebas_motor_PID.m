clear; close all; clc;
Tl = 0.01; % Torque aplicado
torque_values = zeros(1, 6); % Vector para almacenar los valores de torque

for i = 0:5
    Motor_PID_ang;
    Tl = Tl + 0.01;
    torque_values(i+1) = Tl; % Guardar el valor de torque en el vector
end

legend(cellstr(num2str(torque_values', 'Tl = %.2f'))) % Agregar leyenda con los valores de torque
