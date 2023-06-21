% Indice
% x: Estado actual
% k: Paso en el que se encuentra
% accion: Accion de control aplicada
% [y x1 x2]
function [y]=indice_g_motor(x,accion)
   y = x' *diag([1/10 1/25 1/10 44/10])* x+0.001*accion^2;
