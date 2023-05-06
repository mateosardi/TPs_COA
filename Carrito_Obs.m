% TP1_P2_COA
%
% Diseñar un controlador para el péndulo invertido para que evolucione de una posición
% inicial a otra final, analizando desde qué condición inicial de ángulo distinto de cero puede
% controlarse. Considerar que:
% a. Se usa un observador de estados siendo la matriz de salida definida como
% univariable C = [1 0 0 0]

clc;clear all;close all;

% Declaración de constantes
m = 0.1; Fricc = 0.1; long = 0.6; g = 9.8; M = 0.5;
t_etapa = 1e-4; tF = 14;
t=0:t_etapa:tF;
tiempo = round(tF/t_etapa);

% Matrices en el equilibrio inestable
A=[0 1 0 0;0 -Fricc/M -m*g/M 0; 0 0 0 1; 0 Fricc/(long*M) g*(m+M)/(long*M) 0];
B=[0; 1/M; 0; -1/(long*M)];
C=[1 0 0 0];
D=0;
% Controlabilidad y observabilidad
Co = ctrb(A, B);
rank(Co); % = 4 por ende es controlable
Ob = obsv(A,C);
rank(Ob); % = 4 por ende es observable

% Condidiones iniciales
% x = [delta(i); delta_p(i); fi(i); fi_p(i)]
x = [0;0;0;0];
delta=0;
delta_p=0;
fi=0.2;
fi_p=0;
u(1) = 0;
posicion(1) = delta;
angulo(1) = fi;

x_hat=[0;0;0;0]; %Inicializo el Observador

% Diseño de LQR
Q=diag([1 1 10 1]);    R=1e-2;
K = lqr(A,B,Q,R);
eig(A-B*K)

% Controlabilidad y observabilidad
Co = ctrb(A, B);
rank(Co); % = 4 por ende es controlable
Ob = obsv(A,C);
rank(Ob); % = 4 por ende es observable

%Calculo de parametros del observador
Ao=A';
Bo=C';
Co=B';

%Ganancia del observador
Qo=diag([1 1 1 1]); Ro=1e-8;
Ko=(lqr(Ao,Bo,Qo,Ro))';
eig(A-Ko*C)
fi_pp = 0;
delta_pp = 0;

for i=1:1:tiempo+1
    % u=-Ko*x_hat(:,i);  %Con Observador
    u=-K*x;  %Sin Observador
    
    % Ecuaciones diferenciales
    delta_pp = 1/(M+m) *(-m*long*fi_pp*cos(fi)+m*long*(fi_p)^2*sin(fi)-Fricc*delta_p+u);
    fi_pp = (1/long)* (g*sin(fi)-delta_pp*cos(fi));
    
    delta = delta + t_etapa*delta_p;
    delta_p = delta_p+t_etapa*delta_pp;
    fi = fi + t_etapa*fi_p;
    fi_p = fi_p + t_etapa*fi_pp;
    % accion(i) = u;
    % posicion(i) = delta;
    angulo(i) = fi;
    x=[delta; delta_p; fi; fi_p];
end
plot(t,angulo,'r');
hold on;
fi_pp = 0;
delta_pp = 0;
delta=0;
delta_p=0;
fi=0.2;
fi_p=0;
x=[delta; delta_p; fi; fi_p];
for i=1:1:tiempo
    u=-K*x_hat(:,i);  %Con Observador
    y_sal(i)=C*x;
    delta_pp = 1/(M+m) *(-m*long*fi_pp*cos(fi)+m*long*(fi_p)^2*sin(fi)-Fricc*delta_p+u);
    fi_pp = (1/long)* (g*sin(fi)-delta_pp*cos(fi));
    delta = delta + t_etapa*delta_p;
    delta_p = delta_p+t_etapa*delta_pp;
    fi = fi + t_etapa*fi_p;
    fi_p = fi_p + t_etapa*fi_pp;
    x=[delta; delta_p; fi; fi_p];
    %________OBSERVADOR__________
    y_sal_O(i)=C*x_hat(:,i);
    x_hatp = A*x_hat(:,i)+B*u+Ko*(y_sal(i)-y_sal_O(i));
    x_hat(:,i+1)=x_hat(:,i)+t_etapa*x_hatp;
end
plot(t,x_hat(3,:),'b');
