% Indice
% x: Estado actual
% k: Paso en el que se encuentra
% accion: Accion de control aplicada
% [y x1 x2]
function [y]=indice_g(x,accion)

umin=-25; umax=25; pon=0; pon3=0;
% if (accion > umax)
%    pon=5*(accion-umax);
% end
% if (accion < umin)
%    pon=5*abs(umin-accion);
% end


% if x(1) > vfinal
%    pon3= 5*(x(1)-vfinal); %sin sobrepico
% end

% if N==x
%    y = 2000*abs(x(1)-vfinal)+pon+pon3;
% else
   y = x' *diag([1 1 1000 1])* x;
% end
