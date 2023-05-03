close all; clear all; clc

N = 1000;
t = linspace(0,1,N); % Eje temporal

% Señal de pulsos rectangulares
A = 0.5; % Amplitud
T = .1; % Período
d = 20; % Duty cycle

f_t = A/2+A/2*square(2*pi*t/T,d); % Señal f(t)

F_w = fft(f_t);

w = linspace(-pi,pi-2*pi/N,N); % Eje de frecuencia
w = fftshift(w); % Centrar el eje de frecuencia en cero

figure;
plot(t,f_t);
title('f(t)');
xlabel('Tiempo (s)');
ylabel('Amplitud');
figure;
plot(w,F_w);
title('F(W): Espectro de f(t)');
xlabel('Frecuencia (w)');
ylabel('Amplitud');
