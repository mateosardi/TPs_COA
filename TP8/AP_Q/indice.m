% Indice
% x: Estado actual
% k: Paso en el que se encuentra
% accion: Accion de control aplicada
function [y]=indice(k,x,accion)
y=(2+accion)*exp(-x);


