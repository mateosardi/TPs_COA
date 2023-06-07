function [X]=modmotor_PID(t_etapa, xant, accion, Tl)
LAA = 366e-6; J = 5e-9;Ra = 55.6; B = 0; LAF = 5e-6; RF = 64; LFF = 140e-6;
h=1e-7;

ia = xant(1);
i_f = xant(2);
omega = xant(3);
theta= xant(4);

Va = accion;
Vf = 12;

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
