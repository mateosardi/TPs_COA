% Datos obtenidos de:
% International Journal of Dynamics and Control · June 2018
% DOI: 10.1007/s40435-017-0330-x

function [X]=modmotor_PID(t_etapa, xant, accion,Vf, Tl)
LAA = 0.01; J = 0.00208;Ra = 1; B = 0.0011; LAF = 0.1238; RF = 60; LFF = 60;
h=1e-3;

ia = xant(1);
i_f = xant(2);
omega = xant(3);
theta= xant(4);

Va = accion;

for ii=1:t_etapa/h    
    if_p = - RF/LFF * i_f + Vf/LFF;
    ia_p = - Ra/LAA*ia - LAF/LAA*i_f*omega + Va/LAA;
    wp = - B/J*omega + LAF/J * ia * i_f - Tl/J;
    
    ia = ia + h * ia_p;
    i_f = i_f + h * if_p;
    omega = omega + h * wp;
    theta = theta + h * omega;
end
X=[ia, i_f, omega, theta];
end
