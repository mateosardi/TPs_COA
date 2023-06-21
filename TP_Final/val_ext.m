% Funcion que encuentra los valores minimos y maximos de las variables de
% estado
function [min_max] = val_ext(M_max,i_a,i_f,omega,theta)

ia_m = M_max(1,:);
if_m = M_max(2,:);
omega_m = M_max(3,:);
theta_m = M_max(4,:);

if min(i_a)<ia_m(1)
    ia_m(1) = min(i_a);
end
if max(i_a)>ia_m(2)
    ia_m(2) = max(i_a);
end

if min(i_f)<if_m(1)
    if_m(1) = min(i_f);
end
if max(i_f)>if_m(2)
    if_m(2) = max(i_f);
end

if min(omega)<omega_m(1)
    omega_m(1) = min(omega);
end
if max(omega)>omega_m(2)
    omega_m(2) = max(omega);
end

if min(theta)<theta_m(1)
    theta_m(1) = min(theta);
end
if max(theta)>theta_m(2)
    theta_m(2) = max(theta);
end

min_max = [ia_m;if_m;omega_m;theta_m];

end