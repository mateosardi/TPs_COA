% TP5: Dadas las siguientes funciones de transferencia dato, implementar un algoritmo que
% permita obtener sus parámetros estimados a partir de señales de entrada y salida.
% Se pide:
% Emplear como entrada una señal PRBS 7.
% Generar la respuesta al escalón de cada función de transferencia y su correspondiente
% versión identificada.
% Generar la respuesta en frecuencia Bode para cada caso
clc;close all;clear all;

% Parámetros
% orden_a=3; % Orden del denominador
% orden_b=3; % Orden del numerador + 1
% Med=6000;

% Funciones de transferencia
s = tf('s');
F = (2*s-1)/(s^2+3*s+6); ts=1/50; orden_a=2; orden_b=2; Med=6000;
% F = 16*(45*s+1)/((25*s+1)*(30*s+1)); ts=1/20; orden_a=2; orden_b=2; Med=6000;
% F = 1/(0.16*s+1); ts=1/40; orden_a=1; orden_b=1; Med=6000;
% F = (45*s+1)/((0.16*s+1)*(0.4*s^2+0.64*s+1)); ts=1/20; orden_a=3; orden_b=3; Med=6000;

% Entrada escalón
t=0:ts:20000*ts;

StepAmplitude=1.5;
ue = StepAmplitude * square(2*pi*.2*1.0*t);

% Generación de la señal PRBS 7
m=7;
x=ones(m,1);
N=length(ue);%Puntos de la PRBS para muestrear
el=5;
for k=1:el:N
 n_b=xor(x(7),x(6));
 y(k:k+el-1)=x(7);
 for h=m-1:-1:1
 x(h+1)=x(h);
 end
 x(1)=n_b; %Ingreso el nuevo valor
end
x=[];
x=(2*y(1:N)-1);
ue = y(1:N); % Entrada como PRBS

% Obtención de y(t)
[y_D,t_D]=lsim(F,ue,t);
ys=y_D';

%Cuanto es lo que tarda en iniciar el algoritmo
off_set=orden_a+1+50;
u=zeros(Med,1);
z=zeros(Med,1); % es como otra y, auxiliar
u=ue(1+off_set:length(u)+off_set)';ui=u;
z=ys(1+off_set:length(z)+off_set)';zi=z;
for jj=orden_a+1:Med
vec_a=fliplr([ u(jj-orden_b:jj-1); -z(jj-orden_a:jj-1)]');
H(jj-orden_a,:)=vec_a; % Ec 5-8
end
[aa bb ]=size(H);
Z=(z(orden_a+1:end));
in_1=H';
in_2=in_1*H;
in_3=inv(in_2);
in_4=in_3*in_1;
c=in_4*Z; % Ec 5-17, c es el tita que contiene a los a y b

% Evaluación 
abs(roots([1; c(1:orden_a)]))
zo=z;
z=zeros(Med,1);u=ui;
z=zi(1:orden_a);
u=ue(1+off_set:off_set+Med)';
for k=orden_a+1:length(u)
zt=-flip(z(k-orden_a:k-1));
ut=flip(u(k-orden_b:k-1));
z(k)=c'*[zt;ut]; %Ec 125
end
dend=[1; c(1:orden_a)]';numd=[c(orden_a+1:end)]';
sys_id=tf(numd,dend,ts,'Variable','z^-1');
ue=sign(sin(2*pi*.010*t));
[y_sal,t_sal]=lsim(sys_id,ue,t,[0,0]);
[y_D,t_D]=lsim(F,ue,t,[0,0]);
subplot(2,1,1);hold on;
plot(z,'.');
plot(zo,'+r');
plot(ys(1+off_set:off_set+Med),'y');legend('Estimada','Mediciones','Real');legend('Box off'),
title(['Ajuste con el orden b_s/a_s:' num2str(orden_b) '/' num2str(orden_a)]);
xlabel('Muestras')
subplot(2,1,2)
hold on
plot(t_D*ts,y_D,'.');
plot(t_sal*ts,y_sal,'k');legend('Real','Identificada'),legend('Boxoff')
title('Desempeño del modelo ajustado');xlabel('Tiempo. [Seg.]')
figure
step(StepAmplitude*F,'r',StepAmplitude*sys_id,'k'),hold on
legend('Real','Identificada');
legend('boxoff');
F %Sistema original
F_disc = c2d(F,ts);
sysc = d2c(sys_id,'tustin') %sistema identificado
figure
bode(F,logspace(-1,1));
hold on;
bode(sys_id,logspace(-1,1));
legend('Original','Identificada');legend('boxoff')
