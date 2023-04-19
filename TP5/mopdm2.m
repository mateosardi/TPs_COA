%Función carrito
%1998. http://www.math.rutgers.edu/˜sontag/.
function [x] = mopdm2(ts, estado, u)
% x = [delta(i); delta_p(i); fi(i); fi_p(i)]
m=.1;Fricc=0.1; long=0.6;g=9.8;M=.5;
% ts=ts/10;
h=0.0001;
delta = estado(1);
delta_p = estado(2);
fi= estado(3);
fi_p = estado(4);
fi_pp=0;
for i=1:ts/h
    % Ecuaciones diferenciales
    delta_pp = 1/(M+m) *(u-m*long*fi_pp*cos(fi)+m*long*fi_p^2*sin(fi)-Fricc*delta_p);
    fi_pp = (1/long)*(g*sin(fi)-delta_pp*cos(fi));
    delta = delta + h*delta_p;
    delta_p = delta_p+h*delta_pp;
    fi = fi + h*fi_p;
    fi_p = fi_p + h*fi_pp;
end
x=[delta; delta_p; fi; fi_p];