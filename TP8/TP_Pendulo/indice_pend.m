% Indice
% x: Estado actual
% k: Paso en el que se encuentra
% accion: Accion de control aplicada
% Se utiliza el indice de un LQR
function [y]=indice_pend(k,x,accion)
Q=diag([1 1 10 10]);    R=.1;
y = x'*Q*x+accion'*R*accion;
