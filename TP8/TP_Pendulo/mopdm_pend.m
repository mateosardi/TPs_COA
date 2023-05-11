% P?ndulo did?ctico tomado de Sontag Pag 104 Mathematical control theory
%1998. http://www.math.rutgers.edu/?sontag/.
%El controlador es un LQR, s?lo para el punto de equlibrio inestable
%superior, con angulo pr?ximo a cero.
%Autor Juli?n Pucheta.
%?ltima actualizaci?n: 30-05-2020
function [X]=mopdm_pend(tiempo_etapa,xant,accion)
m=.1;Fricc=0.1; long=0.6;g=9.8;M=.5;
tita_pp=0;
h=0.0001;
p=xant(1);
p_p=xant(2);
alfa=xant(3);
omega=xant(4);
for i=1:tiempo_etapa/h
%     estado=[p(i); p_p(i); alfa(i); omega(i)];    
    p_pp=(1/(M+m))*(accion-m*long*tita_pp*cos(alfa)+m*long*omega^2*sin(alfa)-Fricc*p_p);
    tita_pp=(1/long)*(g*sin(alfa)-p_pp*cos(alfa));
    p_p=p_p+h*p_pp;
    p=p+h*p_p;
    omega=omega+h*tita_pp;
    alfa=alfa+h*omega;
end
X=[p; p_p; alfa; omega];
