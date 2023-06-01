% Modelo discreto motor
% k: Tiempo discreto actual
% xant: Estado al que llega
% accion: Accion de control aplicada
% El primer argumento debe ser k.
% El segundo el estado anterior
% El tercero la accion de control
%   3
%-------               Ts=6/109
%s^2+2*s+3
function [X]=mopdm2_motor(tiempo_etapa,xant,accion)
Laa=366e-6; J=5e-9;Ra=55.6;B=0;Ki=6.49e-3;Km=6.53e-3;
Tl =0;
h=1e-7;

ia = xant(1);
w = xant(2);
theta= xant(3);

Va = accion;
ip=0;
for i=1:tiempo_etapa/h
    ip = -Ra/Laa*ia - Km/Laa*theta+Va/Laa;
    w = Ki/J * ia - B/J*theta - Tl/J;
    ia = ia + h*ip;
    theta = theta + h*w;   
end
X=[ia;w;theta];
