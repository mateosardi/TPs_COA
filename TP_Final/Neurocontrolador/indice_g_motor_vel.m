% Indice
% x: Estado actual
% k: Paso en el que se encuentra
% accion: Accion de control aplicada
% [y x1 x2]
function [y]=indice_g_motor(x,accion)
   y = x' *diag([1/1.5 1/1.5 100 1])* x+1.0054*accion^2;
