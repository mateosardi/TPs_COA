% Indice
% x: Estado actual
% k: Paso en el que se encuentra
% accion: Accion de control aplicada
% [y x1 x2]
function [y]=indice_g_motor(x,accion)
   y = x' *diag([0 0 100])* x+0.00000001*accion^2;
