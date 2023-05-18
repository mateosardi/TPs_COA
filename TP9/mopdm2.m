% Modelo discreto pendulo
% k: Tiempo discreto actual
% xant: Estado al que llega
% accion: Accion de control aplicada
% El primer argumento debe ser k.
% El segundo el estado anterior
% El tercero la accion de control
%   3
%-------               Ts=6/109
%s^2+2*s+3
function [X]=mopdm2(tiempo_etapa,xant,accion)
m=.1;Fricc=0.1; long=0.6;g=9.8;M=.5;T0=0.01;
h=0.0001;tiempo=(100/h);p_pp=0;tita_pp=0; t=0:h:tiempo*h;
omega=0:h:tiempo*h; alfa=0:h:tiempo*h; p=0:h:tiempo*h;
p_p=0:h:tiempo*h; u=linspace(0,0,tiempo+1);

p(1)=xant(1);
p_p(1)=xant(2);
alfa(1)=xant(3);
omega(1)=xant(4);
for i=1:tiempo_etapa/h
%     estado=[p(i); p_p(i); alfa(i); omega(i)];
    %     for h_inte=1:/h
    p_pp=(1/(M+m))*(accion-m*long*tita_pp*cos(alfa(i))+m*long*omega(i)^2*sin(alfa(i))-Fricc*p_p(i));
    tita_pp=(1/long)*(g*sin(alfa(i))-p_pp*cos(alfa(i)));
    p_p(i+1)=p_p(i)+h*p_pp;
    p(i+1)=p(i)+h*p_p(i);
    omega(i+1)=omega(i)+h*tita_pp;
    alfa(i+1)=alfa(i)+h*omega(i);
    i=i+1;
end
X=[p(i); p_p(i); alfa(i); omega(i)];
