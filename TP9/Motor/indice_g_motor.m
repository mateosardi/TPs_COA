% Indice
% x: Estado actual
% k: Paso en el que se encuentra
% accion: Accion de control aplicada
% [y x1 x2]
function [y]=indice_g_motor(x,accion)
   y = x' *diag([100 1/900 1/4])* x+1*accion^2;
