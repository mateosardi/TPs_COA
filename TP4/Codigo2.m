close all; clear all; clc;

% Planta
s = tf('s');
Gs = (1+0.16*s)/(1+0.64*s+0.4*s^2);

%Generación de una señal PRBS
m=20;
x=ones(m,1);
N=150000;%Puntos de la PRBS para muestrear
el=5;
for k=1:el:N
 n_b=xor(x(20),x(3));
 y(k:k+el-1)=x(20);
 for h=m-1:-1:1
 x(h+1)=x(h);
 end
 x(1)=n_b; %Ingreso el nuevo valor
end
x=[];
x=(2*y(1:N)-1);

% x=randn(size(x));

%Cálculo de la correlación entre señales digitalizadas
At=1/100;N=length(x);
Tmax=N*At;
t=At:At:Tmax;

% Correlación rxx
rx = xcorr(x,x);

% Transformada de Fourier de la autocorrelación
Sxx = fft(rx); 

% Obtención de y(t)
y_t = lsim(Gs, x, t);

% Correlación ryx
ryx = xcorr(x,y_t);

% Transformada de Fourier de la autocorrelación
Syx = fft(ryx);

% Obtención del modelo F_w
F_w = Syx ./ Sxx;

% Gráficos
M1=fix(0.1*N);%Intervalos de correlacion utiles
fmax=1/(2*At);Af=2*fmax/N;w1=0:Af:fmax-Af;w=2*pi*Af;
Af=2*fmax/M1;
w0=Af:Af:fmax;
w0=w0/(pi);

[mag,phase,wout] = bode(Gs);

figure;
subplot(2,1,1);
semilogx(w0(1:M1/4),20*log10(abs(F_w(1:M1/4))),'.k');
hold on;
bodemag(Gs);
% axis([.1        100         -50          10])
title('Bode aproximado')


subplot(2,1,2);
semilogx(w0(1:M1/4),-180/pi*angle(F_w(1:M1/4)),'.k');
hold on;
plot(wout,squeeze(phase));
axis([.1        100         -135          10])

% figure
% bode(Gs);


% % figure;
% plot(t,x,'.-k');title('Función temporal PRBS 7\Deltat=1.');xlabel('Tiempo [seg.]');
