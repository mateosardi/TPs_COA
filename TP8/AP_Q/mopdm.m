% Modelo discreto
% k: Tiempo discreto actual
% xant: Estado al que llega
% accion: Accion de control aplicada
% El primer argumento debe ser k.
% El segundo el estado anterior
% El tercero la accion de control
function [y]=mopdm(k,xant,accion)
y=xant+(2-2*xant+5/4*xant^2-1/4*xant^3)*accion;